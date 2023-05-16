%%%
% 2023.02.16. (február 16, csütörtök), 10:12
%
% EZT HASZNALD
%
function [T,T_Ms] = Step0_Read_Wastewater_data(xls_name,T,r)
arguments
%%
    xls_name = "Data/NNK_Nyers_uj.xlsx";
    T = [];
    r.Clip = false;
    r.Date_End = datetime("today");
end
%%
%  File: load_szennyviz.m
%  Directory: /home/ppolcz/T/_Epid/Sandbox/2022-07-11-Szennyviz
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. July 11. (2022a)
%

% 2023.02.16. (február 16, csütörtök), 09:52
Varosnevek_angolosan = num2cell({
    'Bekescsaba'      'in B\''ek\''escsaba'
    'Debrecen'        'in Debrecen'
    'Deli'            'in Budapest (South)'
    'Eger'            'in Eger'
    'Eszaki'          'in Budapest (North)'
    'Gyor'            'in Gy\H{o}r'
    'Kaposvar'        'in Kaposv\''ar'
    'Kecskemet'       'in Kecskem\''et'
    'Kozponti'        'in Budapest (city centre)'
    'Miskolc'         'in Miskolc'
    'Nyiregyhaza'     'in Ny\''iregyh\''aza'
    'Pecs'            'in P\''ecs'
    'Salgotarjan'     'in Salg\''otarj\''an'
    'Szeged'          'in Szeged'
    'Szekesfehervar'  'in Sz\''ekesfeh\''erv\''ar'
    'Szekszard'       'in Szeksz\''ard'
    'Szolnok'         'in Szolnok'
    'Szombathely'     'in Szombathely'
    'Tatabanya'       'in Tatab\''anya'
    'Veszprem'        'in Veszpr\''em'
    'Zalaegerszeg'    'in Zalaegerszeg'
    'Pestmegye'       'in Budapest (agglomeration)'
    'Budapestrepter'  'at Budapest Airport (BUD)'    % Nincs lakossag
  % 'Bujak'           'in Buj\''ak'                  % Nincs lakossag
    'Hegyeshalom'     'in Hegyeshalom'               % Nincs lakossag
    'Mosonmagyarovar' 'in Mosonmagyar\''ov\''ar'     % Nincs lakossag
    'Paks'            'in Paks'                      % Nincs lakossag
    },1);
VN_map = containers.Map(Varosnevek_angolosan{:});

% Load xls
warning('off','MATLAB:table:ModifiedAndSavedVarnames')
T_Simitott = readtable(xls_name, ...
    'ReadVariableNames',false, ...
    'Sheet','Simitott', ...
    'Range','2:28', ...
    'VariableNamingRule','modify');

% Nyers adatok betoltese
T_Ms = readtable(xls_name,"Sheet","Meresek");

% 2020. szeptember 3. előtt a méréseket nem tekintem megbízhatónak
% T_Ms(T_Ms.Date <= datetime(2020,09,03),:) = [];

% Torlok minden olyan sort ahol a datum hianyzik
T_Ms(ismissing(T_Ms.Date),:) = [];

% Remove NaNs
T_Ms(isnan(T_Ms.Cp),:) = [];
T_Ms(isnan(T_Ms.Eredmeny_GcpL_szamitott_LOD),:) = [];

% Cq nem lehet nagyobb 40-nel
T_Ms(T_Ms.Cp > 40,:) = [];

% NNK altal kétesnek mínősített PCR eredménynek törlése
T_Ms(contains(T_Ms.Megjegyzes,"kétes PCR eredmény"),:) = [];
T_Ms(contains(T_Ms.Megjegyzes,"alacsony jelintenzitással"),:) = [];

% Nyers adatokban vagas egy adott datumnal (meg mielott barmit is
% elkezdenek csinalni veluk)
if r.Clip
    T_Ms(T_Ms.Date > r.Date_End,:) = [];
end

% Városnevek átnevezése
T_Ms.Varos = categorical(lower(remove_accents_from_string(T_Ms.Varos)));
T_Ms.Szennyviz_Tisztito = categorical(lower(remove_accents_from_string(T_Ms.Szennyviz_Tisztito)));

% Budapest agglomeráció átnevezése
T_Ms.Varos(T_Ms.Varos == "tokol, biatorbagy, szigetszentmiklos, budakeszi, szazhalombatta") = "pestmegye";
T_Ms.Varos(T_Ms.Varos == "budapest repuloter") = "budapestrepter";
% --
ldx = T_Ms.Varos == "budapest";
T_Ms.Varos(ldx) = T_Ms.Szennyviz_Tisztito(ldx);
T_Ms(ismissing(T_Ms.Varos),:) = [];

% Irrelevans varosok adatainak torlese
T_Ms(T_Ms.Varos == "bujak",:) = [];

% Napi szennyviz terfogat hianyzo adatok torlese
% T_Ms(ismissing(T_Ms.Napi_Szennyviz_Terfogat_m3),:) = [];

%% Lakossag kinyerese

Szv_CATS = categories(T_Ms.Varos);

T_Simitott.Var2(1:2) = NaN;
Lakossag = T_Simitott.Var2;
Telep_Idx = find(~isnan(Lakossag));

% (Varos,Lakossag) tabla letrehozasa
Lakossag = Lakossag(Telep_Idx);
Varosok = categorical(lower(remove_accents_from_string(T_Simitott.Var1(Telep_Idx))));
T_Lakossag = table(Varosok,Lakossag,'VariableNames',["Varos","Lakossag"]);

% Budapest agglomeráció átnevezése
T_Lakossag.Varos(T_Lakossag.Varos == "tokol, biatorbagy, szigetszentmiklos, budakeszi, szazhalombatta") = "pestmegye";
T_Lakossag.Varos(T_Lakossag.Varos == "budapest repuloter") = "budapestrepter";

% Nem minden telepulesnek van megadva a lakossaga, de a join miatt ezeket
% is bele kell tenni a tablazatba
ldx = ~ismember(Szv_CATS,T_Lakossag.Varos);
T0_Lakossag = table(categorical(Szv_CATS(ldx)),zeros(sum(ldx),1),'VariableNames',T_Lakossag.Properties.VariableNames);
T_Lakossag = [T_Lakossag ; T0_Lakossag];

% A nyers szennyviz adatokba is bevezetem a relevans lakossag szamokat
T_Ms = join(T_Ms,T_Lakossag);
%{
    % Ha a joint hibat dobna:
    Left = unique(T_Ms(:,"Varos"));
    Right = unique(T_Lakossag(:,"Varos"));
    ism = ismember(Left,Right);
    Left(~ism,:)
%}

% Ahol a lakossag nem adott
%{
    T_Ms_Nincs_lakossag = T_Ms(T_Ms.Lakossag == 0,:);
    unique(T_Ms_Nincs_lakossag.Varos)
%}
% T_Ms(T_Ms.Lakossag == 0,:) = [];

%% Characterize error
% 
% w = T_Ms.Ev*100 + T_Ms.Het;
% 
% [~,Idx] = unique(w);
% 
% Nr = height(T_Ms);
% 
% Idx_From = Idx([1 , 1:end-1]);
% Idx_To = [ Idx(3:end)-1 ; Nr ; Nr ];
% 
% Idx_weekly = cellfun(@(a,b) {a:b},num2cell(Idx_From),num2cell(Idx_To));
% 
% Weekly_wk = cellfun(@(Idx) {T_Ms_nz.Het(Idx)'},Idx_weekly);
% IS_ZERO = norm(cellfun(@(Idx) norm(diff(T_Ms_nz.Het(Idx)')),Idx_weekly));
% 
% Indices, where Cp > 0
% Weekly_ldx = cellfun(@(Idx) {T_Ms.Cp(Idx)' > 0},Idx_weekly);
% Weekly_vals = cellfun(@(Idx,ldx) {T_Ms.Eredmeny_GcpL_szamitott_LOD(Idx(ldx))'},Idx_weekly,Weekly_ldx);
% Weekly_Mean = cellfun(@(v) mean(v,"omitnan"),Weekly_vals);
% Weekly_Std = cellfun(@(v) std(v,"omitnan"),Weekly_vals);
% 
% T_Wkly = table(T_Ms.Ev(Idx),T_Ms.Het(Idx),Weekly_Mean,Weekly_Std,'VariableNames',["Ev","Het","Weekly_Mean","Weekly_Std"]);
% 
% T_Ms = join(T_Ms,T_Wkly);
% 
% LoD = (1/2.5 * 30) * (1000 / 140) * (1000 / 50);
% 
% Error = randn([Nr,1]) .* T_Ms.Weekly_Std;
% T_Ms.Eredmeny_GcpL_szamitott_LOD = max(T_Ms.Eredmeny_GcpL_szamitott_LOD + Error,LoD);
% 
%% Remove outliers as proposed by 2021_Fernandez-Cassi.etal
% Remove records with a higher/lower Cp than the median +- 1.5 cycles 

T_Ms_Cp = table2timetable(unstack(T_Ms(:,["Date","Cp","Varos"]),"Cp","Varos","AggregationFunction",@max));
% T_Ms_Cp_daily = retime(T_Ms_Cp,"daily");
T_Ms_Cp_daily = sortrows(T_Ms_Cp);

isout = isoutlier(T_Ms_Cp_daily,"movmedian",days(30)); % <----- RM OUTLIERS
Vars = T_Ms_Cp_daily.Variables;
Vars(isout == 1) = NaN;
Vars = fillmissing(Vars,"nearest"); % <---- FILL MISSING VALUES
Vars = movmean(Vars,10,1); % <------------- MOOVING MEDIAN FILTER (10 days)
T_Ms_Cp_daily.Variables = Vars;

T_medCp = stack(timetable2table(T_Ms_Cp_daily),T_Ms_Cp_daily.Properties.VariableNames, ...
    'IndexVariableName','Varos','NewDataVariableName','med_Cp');

T_Ms = join(T_Ms,T_medCp,"Keys",["Date","Varos"]);
%{
    % Ha a joint hibat dobna:
    Left = unique(T_Ms(:,["Date","Varos"]));
    Right = unique(T_medCp(:,["Date","Varos"]));
    ism = ismember(Left,Right);
    Left(~ism,:)
%}

% Find outliers, highly deviated measurements from the median Cq
ldx = abs(T_Ms.Cp - T_Ms.med_Cp) > 3;

%%

if ~isempty(T)
    %% Remove outliers
    T_Ms(ldx,:) = [];
else
    %% Remove outliers
    T_Ms_nofilt = T_Ms;
    T_Ms_skipped = T_Ms(ldx,:);
    T_Ms(ldx,:) = [];

    %%

    Tmp = T_Ms_nofilt(:,["Date","Varos","Eredmeny_GcpL_szamitott_LOD"]);

    %% Városonként lineárisan ábrázolni
    
    XLim = [
        datetime(2020,08,31)
        T_Ms.Date(end)
        ];

    FontSize = 12;
    Args = {"Interpreter","latex","FontSize",FontSize};

    % Csak ami kell egy külön táblába:
    R0 = T_Ms_nofilt(:,["Date","Varos","Eredmeny_GcpL_szamitott_LOD","Cp"]);
    R0 = renamevars(R0,"Eredmeny_GcpL_szamitott_LOD","Result");

    Varosok = unique(R0.Varos);

    % Csak ami kell egy külön táblába:
    R1 = T_Ms(:,["Date","Varos","Eredmeny_GcpL_szamitott_LOD","Cp"]);
    R1 = renamevars(R1,"Eredmeny_GcpL_szamitott_LOD","Result");
    
    % Csak ami kell egy külön táblába:
    Tmp = T_Ms_nofilt(:,["Date","Varos","Cp"]);
    T_Ms_Cp_daily_noft = table2timetable(unstack(Tmp,"Cp","Varos","AggregationFunction",@mean));
    
    % Limit of detection
    LoD = (1/2.5 * 30) * (1000 / 140) * (1000 / 50);
    
    Plot_Colors
    
    Col = 7;
    Row = 4;

    if true
        fig = figure(1312);
        fig.Position(3:4) = [3583 1635];
        Tl = tiledlayout(Row*3,Col,'TileSpacing','compact','Padding','compact');
        % Tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
        
        for i = 1:min(numel(Varosok),prod(Tl.GridSize))
            varos = char(Varosok(i));
            Varos = [ upper(varos(1)) varos(2:end) ];

            ldx0 = R0.Varos == varos & XLim(1) <= R0.Date & R0.Date < XLim(2);
            ldx1 = R1.Varos == varos & XLim(1) <= R1.Date & R1.Date < XLim(2);

            ax = nexttile([2 1]); hold on

            plot(R0.Date(ldx0),log10(R0.Result(ldx0)),'.-','MarkerSize',20,'Color',[1 1 1]*0.7)
            plot(R1.Date(ldx1),log10(R1.Result(ldx1)),'.-','MarkerSize',20,'Color',Color_1)
            title("PCR results " + VN_map(Varos),Args{:})

            plot(XLim,[1 1]*log10(LoD),'Color',Color_5,'LineWidth',2)

            grid on, box on
            xlim(XLim)

            if mod(i,Col) == 1
                ylabel('log GC/L')
            end

            drawnow

            Idx = tilenum(ax) + 2*Col;
            ax = nexttile(Idx); hold on

            Cp = T_Ms_Cp_daily_noft.(varos);

            ldx = ~isnan(Cp);

            plot(R0.Date(ldx0),R0.Cp(ldx0),'.-','MarkerSize',20,'Color',[1 1 1]*0.7)
            plot(R1.Date(ldx1),R1.Cp(ldx1),'.-','MarkerSize',20,'Color',Color_1)
            plot(T_Ms_Cp_daily.Date,T_Ms_Cp_daily.(varos),'LineWidth',3,'Color',Color_2)
            title("$C_\mathrm{q}$ " + VN_map(Varos),Args{:})

            grid on, box on
            xlim(XLim)
            ylim([0 40])
            if mod(i,Col) == 1
                ylabel('cycle (Cq)')
            end

            drawnow
        end
    end

    DIR = "/home/ppolcz/T/Dropbox/Peti/Munka/01_PPKE_2020/Dokumentaciok/Docs_CsutakB_PhD/07_COVID-Szennyviz/actual/fig";
    clipboard("copy",DIR)
    fname = DIR + "/PCR-Results_All.pdf";
    % exportgraphics(fig,fname,"ContentType","vector")
    
    %%

    XLim = [
        datetime(2020,08,31)
        T_Ms.Date(end)
        ];
    Date = XLim(1):XLim(2);
    
    FontSize = 12;
    Args = {"Interpreter","latex","FontSize",FontSize};


    if true
        
        for i = 1:numel(Varosok)
            varos = char(Varosok(i));
            %%
            % varos = 'szeged';
            Varos = [ upper(varos(1)) varos(2:end) ];

            fig = figure(734);
            fig.Position(3:4) = [619 447];
            Tl = tiledlayout(3,1,"TileSpacing","tight","Padding","compact");

            ldx0 = R0.Varos == varos & XLim(1) <= R0.Date & R0.Date < XLim(2);
            ldx1 = R1.Varos == varos & XLim(1) <= R1.Date & R1.Date < XLim(2);

            Ax = nexttile([2 1]); hold on, grid on
            Ax.YScale = 'log';
            % Xl = xline(Date(day(Date) == 1),'Color',[1 1 1]*0.75);
            % Ax.YGrid = "on";
            % for xl = Xl'
            %     xl.HandleVisibility = "off";
            % end

            plot(R0.Date(ldx0),R0.Result(ldx0),'.-','MarkerSize',20,'Color',[1 1 1]*0.7,...
                'DisplayName','Outliers')
            plot(R1.Date(ldx1),R1.Result(ldx1),'.-','MarkerSize',20,'Color',Color_1,...
                'DisplayName','Accepted measurments')
            % title("PCR results " + VN_map(Varos) + " \rule[1em]{0.001pt}{0.001pt}",Args{:});
            
            plot(XLim,[1 1]*LoD,'Color',Color_5,'LineWidth',2,...
                'DisplayName','Limit of detection')

            box on
            xlim(XLim)
            ylim([10e2 10e6])
            ylabel('Genome copy concentration [GC/L]',Args{:})
            Leg = legend('Location','northwest',Args{:});

            Ax.XTickLabel = [];

            drawnow

            Ax(2) = nexttile; hold on, grid on
            % Xl = xline(Date(day(Date) == 1),'Color',[1 1 1]*0.75);
            % Ax(2).YGrid = "on";
            % for xl = Xl'
            %     xl.HandleVisibility = "off";
            % end

            Cp = T_Ms_Cp_daily_noft.(varos);

            ldx = ~isnan(Cp);

            plot(R0.Date(ldx0),R0.Cp(ldx0),'.-','MarkerSize',20,'Color',[1 1 1]*0.7,'HandleVisibility','off')
            plot(R1.Date(ldx1),R1.Cp(ldx1),'.-','MarkerSize',20,'Color',Color_1,'HandleVisibility','off')
            plot(T_Ms_Cp_daily.Date,T_Ms_Cp_daily.(varos),'LineWidth',3,'Color',Color_2,'DisplayName','Moving median')
            Leg = legend('Location','east',Args{:});

            box on
            xlim(XLim)
            ylim([0 40])
            ylabel('Cycles ($C_\mathrm{q}$)',Args{:})

            drawnow

            for idx = 1:numel(Ax)
                ax = Ax(idx);
                Logger.latexify_axis(ax,FontSize);
                % ax.XTick = Date(day(Date) == 1 & mod(month(Date),3) == 1);
                ax.XTick = Date(day(Date) == 1 & mod(month(Date),6) == 3);
                ax.XMinorGrid = 'on';
                ax.XAxis.MinorTick = 'off';
                ax.XAxis.MinorTickValues = Date(day(Date) == 1);

                drawnow
            end

            DIR = "/home/ppolcz/T/Dropbox/Peti/Munka/01_PPKE_2020/Dokumentaciok/Docs_CsutakB_PhD/07_COVID-Szennyviz/actual/fig";
            clipboard("copy",DIR)
            tag = '';
            if isempty(ax.Title.String)
                tag = 'notit-';
            end
            fname = DIR + "/PCR-Results_" + Varos + "-" + tag + string(datetime(date,'Format','uuuu-MM-dd')) + ".pdf";
            % exportgraphics(fig,fname,"ContentType","vector")

        end
    end

    %% NNK-s módszerrel számolom ki a szennyvíz görbét
    
    % Csak ami kell egy külön táblába:
    Tmp = T_Ms(:,["Date","Ev","Het","Varos","Eredmeny_GcpL_szamitott_LOD","Lakossag"]);
    Tmp = renamevars(Tmp,"Eredmeny_GcpL_szamitott_LOD","Result");
    Tmp.LogResult = log10(Tmp.Result);
    Tmp.WghtLogResult = Tmp.LogResult .* Tmp.Lakossag;
    
    % Súlyozott átlag kiszámítása
    S = groupsummary(Tmp,["Ev","Het","Varos"],["sum","mean"],["Result","Lakossag","LogResult","WghtLogResult"]);
    S_Date = groupsummary(Tmp,["Ev","Het","Varos"],["min"],"Date");
    S = join(S_Date,S);
    S.min_Date.Format = 'uuuu-MM-dd (eee)';
    
    % Átlag dátumát a hét keddi napjára igazítom
    S.Date = S.min_Date - weekday(S.min_Date)+3;
    S.Het = week(S.Date);
    S.Date.Format = 'uuuu-MM-dd (eee)';
    S = movevars(S,"Date","Before","Ev");
    
    % Az átlagoló táblából kiszedem azokat, amelyekkel lehet unstack
    % segítségével városonként oszlopokba rendezni
    S1 = S(:,["Date","Varos","mean_LogResult"]);
    S1.Date.Format = "'D'_uuuu_MM_dd";
    
    % Megtörténik a váronoskénti adatrendezés
    US = unstack(S1,"mean_LogResult","Date");
    
    % Ismét bevezetem a lakosságot a második oszlopba
    US = join(US,T_Lakossag);
    US = movevars(US,"Lakossag","After","Varos");
    
    Date = datetime(US.Properties.VariableNames(3:end)',"InputFormat","'D'_uuuu_MM_dd");
    
    % Elso sorba bevezetem a hetek számát (miheztartás végett, illetve, hogy
    % nagyon hasonlítson az NNK-s táblához)
    US = US([1 1:end],:);
    US.Varos(1) = "Week";
    US.Lakossag(1) = NaN;
    Wt = US(1,3:end);
    Wt.Variables = week(Date)';
    US(1,3:end) = Wt;
    
    % Háromheti mozgótátlagot számolok városonként
    Wt = US(2:end,3:end);
    Wt.Variables = movmean(Wt.Variables,3,2,'omitnan');
    US_ma = US;
    US_ma(2:end,3:end) = Wt;
    
    % Súlyozott átlagot számolok a városok simított adatsorából
    WghtLogAverage = sum(US.Lakossag(2:end) .* Wt.Variables,'omitnan') ./ sum(US.Lakossag(2:end) .* ~isnan(Wt.Variables));
    Result = 10.^WghtLogAverage';
    
    % Vegül elmentem egy táblázatba
    T_NNKrek = timetable(Date,Result);
    clear Date Result
end

%%

T_Simitott.Var2 = [];

% These values are written in logarithms
LogData = table2array(T_Simitott(Telep_Idx,2:end));

%%%
% 2023.01.07. (január  7, szombat), 04:54
    
    Year = table2array(T_Simitott(1,2:end))';
    Week = table2array(T_Simitott(2,2:end))';
    
    % Fill missing years
    if isnan(Year(1))
        Year(1) = 2020;
    end
    for i = 2:numel(Year)
        if isnan(Year(i))
            Year(i) = Year(i-1);
        end
    end
    
    Date_Monday = pcz_week2date(Year,Week);
    Date_Monday.Format = "uuuu-MM-dd (eeee)";

%%%

% Linear-scale values
% Data = 10.^LogData;
Data = LogData;

% The weighted average is computed in logirithmic scale
Relevans_Lakossag = ~isnan(Data) .* Lakossag;
Data(isnan(Data)) = 0;
Result = sum(Data .* Relevans_Lakossag,1)' ./ sum(Relevans_Lakossag,1)';

% Linear-scale values
Result = 10.^Result;

% 2023.01.06. (január  6, péntek), 13:25
Date = Date_Monday + 1; % Tuesdays

O = timetable(Date,Date_Monday,Year,Week,Result);

O.Linear = Data';
O.Log = log10(Data');

N = days(O.Date(end) - O.Date(1));
Days = (0:N)';
Days_raw = days(O.Date - O.Date(1));
O_Date = Days + O.Date(1);
ww = spline(Days_raw,O.Result,Days);
wwLin = interp1(Days_raw,O.Result,Days);

%%
%%%

% Csak ami kell egy külön táblába:
Tmp = T_Ms(:,["Date","Eredmeny_GcpL_szamitott_LOD","Lakossag"]);
Tmp = renamevars(Tmp,"Eredmeny_GcpL_szamitott_LOD","Result");

% Sulyozott atlag nevezoje: napokra osszlakossag
Date = Tmp.Date;
Szv = timetable(Date);
Szv.Lakossag = Tmp.Lakossag;

% Sulyozott atlag szamlaloja: eredmeny*lakossag osszege
Szv.WghtSum = Tmp.Result .* Tmp.Lakossag;
Szv = retime(Szv,"daily","sum");

% Kiszamolom a sulyozott atlagot
Szv(Szv.Lakossag == 0,:) = [];
Szv.Value = Szv.WghtSum ./ Szv.Lakossag;
Szv.Value = medfilt1(Szv.Value,5);

Szv(:,["Lakossag","WghtSum"]) = [];

% Megjegyzem, hogy melyik adatpont volt az utolso
Date_Last_Ms = Szv.Date(end);

% Linearisan interpolalom a hianyzo napokat
Szv = retime(Szv,"daily","spline");

% Egy kis csalas: az utolso het adatait megismetlem ketszer, hogy a szuro
% eredmenye ne konvergaljon nullaba, hanem tartsa az utolso het szamait.
for i = 1:2
    Szv = tb_expand_timerange(Szv,[Szv.Date(1) Szv.Date(end)+7]);
    Szv(end-6:end,:) = Szv(end-13:end-7,:);
end

%%

CutOff = 1/21;
CO_Perc = CutOff / (1/2);
[b,a] = butter(4,CO_Perc);
Szv.Filt = filtfilt(b,a,Szv.Value);

% Az utolso ket hetet, amelyet csak a szures miatt toldottam oda, most
% kiveszem.
Szv(Szv.Date > Date_Last_Ms,:) = [];

if isempty(T) && false
%%
    XLim = [
        datetime(2020,08,01)
        Szv.Date(end)
        ];
    
    fig = figure(153);
    fig.Position(3:4) = [595 303];
    Tl = tiledlayout(1,1,"TileSpacing","compact","Padding","compact");
    
    % N = 1024*4;
    % F = abs(fft(Szv.Value,N));
    % ax = nexttile; hold on; grid on; box on
    % plot(linspace(0,0.5,N/2),F(1:N/2));
    % FTicks = sort(1./[365 61 30 21 14 7:-1:2 , 3.5 , 2.33]);
    % FTickLabels = cellfun(@(n) {num2str(n)},num2cell(1./FTicks));
    % ax.XTick = FTicks;
    % ax.XTickLabel = FTickLabels;
    % xlabel('Signal period in days')
    % xline(CutOff,'r-',sprintf('Cut-off frequency (%d days: %2.f%%)',1/CutOff,CO_Perc*100),'LabelOrientation','horizontal')
    
    ax = nexttile; hold on; grid on; box on
    title 'Filtered in linear scale'

    if exist('T_Data','var')
        yyaxis right
        plot(T_Data.Date,T_Data.H_off_ma,'k-','LineWidth',2, ...
        'DisplayName','K\''orh\''azban \''apolt (hivatalos adatsor)')
        ylim([0 15000])
        ylabel('K\''orh\''azban \''apolt','Interpreter','latex')

        yyaxis left
    end

    Pl = plot(T_Ms.Date,T_Ms.Eredmeny_GcpL_szamitott_LOD,'.','Color',[1 1 1]*0.8,'MarkerSize',15, ...
        'DisplayName','Raw measurements (filtered according to Fernandez-Cassi et.al. (2021) from Switzerland');
    Pl = plot(Szv.Date,Szv.Value,'-','Color',Color_1, ...
        'DisplayName','Averaged and interpolated series');
    Pl = plot(Szv.Date,Szv.Filt,'-','Color',Color_2,'LineWidth',2, ...
        'DisplayName','Zero-phase filtered series ($C_k^\mathrm{Off}$)');
    Pl = plot(O_Date,wwLin,'-','Color',Color_4,'LineWidth',2, ...
        'DisplayName','A sz\H{u}retlen nyers adat \''atlagol\''asa \''es sim\''it\''asa logaritmus tartom\''anyban (NNK)');
    plot(T_NNKrek.Date,T_NNKrek.Result,'-','Color',Color_5,'LineWidth',2, ...
        'DisplayName','A SZ\H{U}RT nyers adat \''atlagol\''asa \''es sim\''it\''asa logaritmus tartom\''anyban')
    
    Font = {'Interpreter','latex','FontSize',12};
    
    Logger.latexify_axis(ax,12)
    ylabel('Average genome copy conc. [GC/L]',Font{:})
    ax.YLim(1) = 0;

    xlim(XLim)
    ylim([0,1.5e6])

    Date = datetime(2021,10,15,'Format','uuuu-MM-dd');
    Xl = xline(Date,'-',"Delta (" + string(Date) + ")",'Interpreter','latex','FontSize',14,'HandleVisibility','off');

    Date = datetime(2022,01,01,'Format','uuuu-MM-dd');
    xline(Date,'-',"Omicron BA.1 (" + string(Date) + ")",'Interpreter','latex','FontSize',14,'HandleVisibility','off')

    Date = datetime(2022,03,01,'Format','uuuu-MM-dd');
    xline(Date,'-',"Omicron BA.2 (" + string(Date) + ")",'Interpreter','latex','FontSize',14,'HandleVisibility','off')

    Leg = legend(Font{:},'Location','northwest','Box','off');

    % DIR = "/home/ppolcz/T/Dropbox/Peti/Munka/01_PPKE_2020/Dokumentaciok/Docs_CsutakB_PhD/07_COVID-Szennyviz/actual/fig";
    % clipboard("copy",DIR)
    % fname = DIR + "/Raw_vs_filtered_wwdata.pdf";
    % exportgraphics(ax,fname,"ContentType","vector")

    return
end

if isempty(T)
%%
    Plot_Colors

    XLim = [
        datetime(2021,09,30)
        datetime(2022,10,01)
        ];

    fig = figure(153);
    fig.Position(3:4) = [595 303];
    Tl = tiledlayout(1,1,"TileSpacing","compact","Padding","compact");

    % N = 1024*4;
    % F = abs(fft(Szv.Value,N));
    % ax = nexttile; hold on; grid on; box on
    % plot(linspace(0,0.5,N/2),F(1:N/2));
    % FTicks = sort(1./[365 61 30 21 14 7:-1:2 , 3.5 , 2.33]);
    % FTickLabels = cellfun(@(n) {num2str(n)},num2cell(1./FTicks));
    % ax.XTick = FTicks;
    % ax.XTickLabel = FTickLabels;
    % xlabel('Signal period in days')
    % xline(CutOff,'r-',sprintf('Cut-off frequency (%d days: %2.f%%)',1/CutOff,CO_Perc*100),'LabelOrientation','horizontal')

    ax = nexttile; hold on; grid on; box on
    % title 'Filtered in linear scale'

    if exist('T_Data','var')
        yyaxis right
        plot(T_Data.Date,T_Data.H_off_ma,'k-','LineWidth',2, ...
        'DisplayName','K\''orh\''azban \''apolt (hivatalos adatsor)')
        ylim([0 15000])
        ylabel('K\''orh\''azban \''apolt','Interpreter','latex')

        yyaxis left
    end

%     Pl = plot(T_Ms_skipped.Date,T_Ms_skipped.Eredmeny_GcpL_szamitott_LOD,'.','Color',[1 1 1]*0.8,'MarkerSize',15, ...
%         'DisplayName','Raw measurements ($C_\mathrm{raw}$)');
    Pl = plot(T_Ms.Date,T_Ms.Eredmeny_GcpL_szamitott_LOD,'.','Color',[1 1 1]*0.8,'MarkerSize',15, ...
        'DisplayName','Non--inhibited measurements ($C_\mathrm{raw}$)');
    Pl = plot(Szv.Date,Szv.Value,'-','Color',Color_1, ...
        'DisplayName','Averaged and interpolated series');
    Pl = plot(Szv.Date,Szv.Filt,'-','Color',Color_2,'LineWidth',2, ...
        'DisplayName','Zero-phase filtered series ($C_k^\mathrm{Off}$)');
    Pl = plot(O_Date,wwLin,'-','Color',Color_4,'LineWidth',2, ...
        'DisplayName','Filtered by R\''oka et al. (2021)');
%     plot(T_NNKrek.Date,T_NNKrek.Result,'-','Color',Color_5,'LineWidth',2, ...
%         'DisplayName','A SZ\H{U}RT nyers adat \''atlagol\''asa \''es sim\''it\''asa logaritmus tartom\''anyban')

    Font = {'Interpreter','latex','FontSize',12};

    Logger.latexify_axis(ax,12)
    ylabel('Average genome copy conc. [GC/L]',Font{:})
    ax.YLim(1) = 0;

    xlim(XLim)
    ylim([0,1.5e6])

    Leg = legend(Font{:},'Location','northeast','Box','off');

    % %{
        DIR = "/home/ppolcz/T/Dropbox/Peti/Munka/01_PPKE_2020/Dokumentaciok/Docs_CsutakB_PhD/07_COVID-Szennyviz/actual/fig";
        clipboard("copy",DIR)
        fname = DIR + "/Raw_vs_filtered_wwdata.pdf";
        exportgraphics(ax,fname,"ContentType","vector")
    %}
    return
end

%% Bevezetem a nagy tablazatba

T = tb_expand_timerange(T,Szv.Date([1,end]));

Idx = Szv.Date(1) <= T.Date & T.Date <= Szv.Date(end);
T.Szennyviz(Idx) = max(Szv.Filt,0);
T.Szennyviz_Nyers(Idx) = Szv.Value;
T.Properties.UserData.Date_Last_Szennyviz = Date_Last_Ms;

end
