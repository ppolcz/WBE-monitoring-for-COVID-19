function T = epid_load_cowid(XLS_NAME,r)
arguments
    XLS_NAME;
    r.Clip = false;
    r.Date_End = datetime("today");
end
%%
%  File: pf_load_data.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study7_tracking
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2021. March 31. (2020b)
%  Major revision on 2022. May 05. (2022a)
%  Major revision on 2023. April 19. (2022b)
%

try
    %% Download and load csv from Our World in Data

    opts = detectImportOptions(XLS_NAME);
    opts = setvartype(opts,"iso_code","categorical");
    opts = setvartype(opts,"date","datetime");
    opts = setvartype(opts,[ ...
        "hosp_patients" ...
        "new_cases" ...
        "total_deaths" ...
        "people_vaccinated" ...
        "people_fully_vaccinated" ...
        "total_boosters" ...
        "stringency_index" ...
        "new_tests" ...
        "positive_rate"],"double");
    
    Main = readtable(XLS_NAME,opts,'ReadVariableNames',true);

    % 2022.12.15. (december 15, csütörtök), 15:31
    if r.Clip
        Main(Main.date > r.Date_End,:) = [];
    end

    %% Remove textual variables

    Main(:,vartype('string')) = [];
    Main(:,vartype('cellstr')) = [];
    Main(:,vartype('categorical')) = [];

    %% Convert it to a timetable

    Main = table2timetable(Main,"RowTimes","date");
    
    %% Interpolate data

    if isnan(Main.Properties.TimeStep) ...
        || (isa(Main.Properties.TimeStep,"duration") && days(Main.Properties.TimeStep) ~= 1) ...
        || (isa(Main.Properties.TimeStep,"calendarDuration") && caldays(Main.Properties.TimeStep) ~= 1)

        Main = retime(Main,'daily','spline');
    end

    %%
    
    % Load vaccination information
    try
        V1 = resolve_missing_weekends(Main.people_vaccinated);
    catch
        V1 = zeros(height(Main),1);
    end

    try
        V2 = resolve_missing_weekends(Main.people_fully_vaccinated);
    catch
        V2 = zeros(height(Main),1);
    end

    if sum(~isnan(Main.total_boosters)) < 2
        Vb = zeros(size(Main.date));
    else
        Vb = resolve_missing_weekends(Main.total_boosters);
    end
    
    
    % Load stringency index
    StrIdx = (100 - resolve_missing_weekends(Main.stringency_index))/100;
    
    %% Date variables
    
    Start_Date = Main.date(1);
    End_Date = Main.date(end);
    
    %% Simulation time span
    
    N = days(End_Date - Start_Date);
    t = (0:N)';
    assert(size(Main,1) == N+1,'Polcz: Days are missing in OWID''s data table!');
    
    %% Construct data table
    % 
    % In the table the `Date` is the absolute time label.
    % The `Day` or `Day_MPC`, etc. are all relative time labels.
    
    % Window radios and length
    wr = 3;
    w = wr*2+1;
    
    T = table;
    T.Date = Main.date;     % Absolute time label
    T.Day = t;              % Relative time label
    T.Np = Main.population;
    T.D_off = Main.total_deaths;
    T.H_off = Main.hosp_patients;
    T.H_off_ma = movmean( resolve_missing_weekends(Main.hosp_patients,Fill=0), w);
    T.New_Tests = resolve_cumulative_weekends(Main.new_tests);
    T.New_Cases = movmean( resolve_cumulative_weekends(Main.new_cases), w);
    T.Positive_Rate = Main.positive_rate;
    T.V_first = V1;
    T.V_second = V2;
    T.V_boosted = Vb;
    T.Stringency = StrIdx;
    T.Rt = Main.reproduction_rate;

    T = table2timetable(T,"RowTimes","Date");

    % Fill missing values if relevant
    T(:,"D_off") = fillmissing(T(:,"D_off"),'constant',0);
    T(:,"Rt") = fillmissing(T(:,"Rt"),"nearest");
    
    %% Detect last available hospitalization data
    
    Idx = find(~isnan(T.H_off),1,"last");
    T.H_off_ma(Idx+wr+1:end) = NaN;
    
    Idx = min(Idx+wr,size(T,1));
    T.Properties.UserData.Date_Last_Available_H = T.Date(Idx);
    T.Properties.UserData.Date_Last_OWID = End_Date;
    
    %%
    
    % Erre azert van szukseg, hogy ha majd ki kell egisziteni a tablazatot a
    % jovo iranyaba, akkor az utolso ertek, amelyet ismetelgetni fog, az
    % megfelelo legyen. A New_Cases eseten ez 'NaN' kell legyen.
    T = tb_expand_timerange(T,[T.Date(1) T.Date(end)+1]);
    T.New_Cases(end) = NaN;

catch e
    getReport(e)
    T = [];
end

%% You may store further information in `UserData`
% T.Properties.UserData;

end

function X = resolve_cumulative_weekends(X)
%%

    X_old = X;

    X(X == 0) = NaN;

    if ~any(isnan(X))
        return
    end

    we = isnan(X);
    cumulative = [false ; we(2:end) < we(1:end-1)];

    cswe = cumsum(we);
    cswe_cum = cswe(cumulative);

    Idx = find(cumulative);
    X_we = X(cumulative);
    we_length = [cswe_cum(1) ; diff(cswe_cum)];

    % [ we , cswe , cumulative , X ]

    for i = 1:numel(Idx)
        X(Idx(i)-we_length(i):Idx(i)) = X_we(i) / (we_length(i)+1);
    end

    Idx = find(~isnan(X),1,'last');
    X(Idx+1:end) = 0;

    assert(sum(X_old(~isnan(X_old))) - sum(X) < 1e-8, ...
        'Numerical errors in function resolve_cumulative_weekends.');

end

function X = resolve_missing_weekends(X,args)
arguments
    X
    args.Fill = 0;
    args.FillEnd = 'repeat';
    args.Method = 'spline';
end
%%
    
    t = 1:numel(X);

    X_isnum = ~isnan(X);

    t_ref = t(X_isnum);
    X_ref = X(X_isnum);
    
    Idx1 = find(X_isnum,1,'first');
    Idx2 = find(X_isnum,1,'last');

    t_ = Idx1:Idx2;
    X(t_) = interp1(t_ref,X_ref,t_,args.Method);

    if strcmp(args.FillEnd,'repeat')
        X(Idx2+1:end) = X(Idx2);
    end

    X(isnan(X)) = args.Fill;

end
