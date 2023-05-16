%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2023. April 25. (2022b)
% 

Plot_Colors

xls_name = "Data/NNK_Nyers_uj.xlsx";
T = [];
r.Clip = false;
r.Date_End = datetime("today");

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
T_Ms = renamevars(T_Ms,["Koncentralt_Mennyiseg_ml","Koncentratum_Terfogat_ml","Minta_NS_Kivonashoz","Elualas_TF"],["V0","Vc","Vs","Ve"]);

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

%%

hyp_Cp_pos_file = "Results/hyp_Cp_pos.mat";
hyp_Cp_Res_file = "Results/hyp_Cp_Res.mat";

Start_Date = min(T_Ms.Date);
T_Ms.Day = days(T_Ms.Date - Start_Date);

ldx = T_Ms.Cp > 0 & log10(T_Ms.PCR_Eredmeny) >= 0;
if ~exist(hyp_Cp_Res_file,'file')

    hyp_init = struct('mean',[],'cov',zeros(1,2),'lik',-1);
    hyp_Cp_Res = gpml_minimize(hyp_init,@gp,-200,@infGaussLik,[],@covSEard,@likGauss,T_Ms.Cp(ldx),log10(T_Ms.PCR_Eredmeny(ldx)));

    save(hyp_Cp_Res_file,"hyp_Cp_Res")
else
    load(hyp_Cp_Res_file)
end

Cps = linspace(26,40,100);
hyp_Cp_Res = GP_hyp(hyp_Cp_Res,T_Ms.Cp(ldx),log10(T_Ms.PCR_Eredmeny(ldx)));
GP_eval(hyp_Cp_Res);
[GP_mean,GP_var] = GP_eval(hyp_Cp_Res,Cps);

Cq_Opt = 37.2;
XData = 25:0.1:45;
yData = -XData*log10(2) + log10(2^Cq_Opt);
hyp_Cp_Res.X = XData';
hyp_Cp_Res.y = yData';

GP_eval(hyp_Cp_Res);
[GP_mean_,GP_var_] = GP_eval(hyp_Cp_Res,Cps);

%{
%%
fig = figure(10); 
delete(fig.Children);
ax = axes(fig);
hold on; grid on; box on
PlMsm = plot(T_Ms.Cp(ldx),log10(T_Ms.PCR_Eredmeny(ldx)),'.');
Sh = plot_mean_var(Cps,GP_mean_,sqrt(GP_var_),Color_2,"LineWidth",3);
ylim([0,4])
xlim([28,38])
% plot(ax.XLim,-ax.XLim*log10(2) + log10(2^Cq_Opt))
title('$C_\mathrm{q}$ vs. log10(PCR result)','Interpreter','latex','FontSize',14)
xlabel('$C_\mathrm{q}$','Interpreter','latex','FontSize',14)
ylabel('log10(PCR results)','Interpreter','latex','FontSize',14)
legend([PlMsm Sh(1) Sh(4)],{
    'Measurements'
    'Mean'
    '0.95\% CI'
    },"Location","northeast",'Interpreter','latex','FontSize',14)
Logger.latexify_axis(ax,14);
exportgraphics(ax,"/home/ppolcz/Dropbox/Apps/Overleaf/COVID-22-Szennyviz-Revision1-Proba1/fig/GP_Cq_vs_log10PCR_result.pdf","ContentType","vector")
%%
%}

%{
%%
fig = figure(71); 
delete(fig.Children);
ax = axes(fig);
hold on; grid on; box on
PlMsm = plot(log10(T_Ms.PCR_Eredmeny(ldx)),T_Ms.Cp(ldx),'.');
PlMean = plot(GP_mean_,Cps,"LineWidth",3);
xlim([0,3])
ylim([30,40])
Pl35 = plot(-ax.YLim*log10(2) + log10(2^35),ax.YLim,"LineWidth",3);
Pl40 = plot(-ax.YLim*log10(2) + log10(2^40),ax.YLim,"LineWidth",3);
title('log10(PCR result) vs. $C_\mathrm{p}$','Interpreter','latex','FontSize',14)
xlabel('log10(PCR results)','Interpreter','latex','FontSize',14)
ylabel('$C_\mathrm{p}$','Interpreter','latex','FontSize',14)
legend([PlMsm Pl40 PlMean Pl35],{
    'Measurements'
    'Upper: $40\!-\!3.3219\,x$'
    'Opt. fit: $37.2\!-\!3.3219\,x$'
    'Lower: $35\!-\!3.3219\,x$'
    },"Location","northeast",'Interpreter','latex','FontSize',12)
Logger.latexify_axis(ax,14);
exportgraphics(ax,"/home/ppolcz/Dropbox/Apps/Overleaf/COVID-22-Szennyviz-Revision1-Proba1/fig/GP_log10PCR_result_vs_Cq.pdf","ContentType","vector")
%%
%}

%%

if ~exist(hyp_Cp_pos_file,'file')
    hyp_init = struct('mean',[],'cov',zeros(1,2),'lik',-1);
    hyp_Cp_pos = gpml_minimize(hyp_init,@gp,-200,@infGaussLik,[],@covSEard,@likGauss,T_Ms.Day(ldx),T_Ms.Cp(ldx));
    save(hyp_Cp_pos_file,"hyp_Cp_pos")
else
    load(hyp_Cp_pos_file)
end

Days = 0:max(T_Ms.Day);
hyp_Cp_pos = GP_hyp(hyp_Cp_pos,T_Ms.Day(ldx),T_Ms.Cp(ldx));
GP_eval(hyp_Cp_pos);
[Cp_mean,Cp_var] = GP_eval(hyp_Cp_pos,Days);

T_Ms.Pos = double(T_Ms.Cp > 0);
T_Ms.Pr_Pos = zeros(height(T_Ms),1);

r = 7;
T_Ms.Covariance = zeros(height(T_Ms),3);
for i = 1:height(T_Ms)
    ldx = T_Ms.Date(i) - r <= T_Ms.Date & T_Ms.Date <= T_Ms.Date(i) + r;
    T_Ms.Pr_Pos(i) = mean(T_Ms.Pos(ldx));
end
clear r

r = 30;
T_Ms.Covariance = zeros(height(T_Ms),3);
for i = 1:height(T_Ms)
    ldx = T_Ms.Date(i) - r <= T_Ms.Date & T_Ms.Date <= T_Ms.Date(i) + r;
    ldx = ldx & T_Ms.Cp > 0;
    if sum(ldx) > 2
        Sigma = cov([T_Ms.Cp(ldx),log10(T_Ms.PCR_Eredmeny(ldx))]);
        T_Ms.Covariance(i,:) = Sigma([1 2 4]);
    end
end
clear r

%{
%%
fig = figure(11); 
fig.Position(3:4) = [1375 426];
delete(fig.Children);
ax = axes(fig);
hold on;
PlMsm = plot(T_Ms.Date,T_Ms.Cp,'o','MarkerSize',5);
Sh = plot_mean_var(Start_Date+Days,Cp_mean,sqrt(Cp_var),Color_2,"LineWidth",3);
grid on; box on
ax.XLim = Start_Date + [0 Days(end)];
ax.YLim = [-0.1 40];
Logger.latexify_axis(ax,14);
ylabel('$C_\mathrm{q}$ ($0$ for negative tests)','Interpreter','latex','FontSize',16)
title('$C_\mathrm{q}$ values and the probability of a test being positive','Interpreter','latex','FontSize',16)

[~,Idx] = unique(T_Ms.Date);

XData_ = T_Ms.Date(Idx);
YData_ = T_Ms.Pr_Pos(Idx)*100;

XData = XData_(1):XData_(end);
YData = interp1(XData_,YData_,XData);

yyaxis right
PlPr = bar(XData,YData,1,'FaceColor',Color_6,'FaceAlpha',0.5);
plot(T_Ms.Date([1,end]),[1 1]*0,'--','Color',[1 1 1]*0.5)
plot(T_Ms.Date([1,end]),[1 1]*50,'--','Color',[1 1 1]*0.5)
plot(T_Ms.Date([1,end]),[1 1]*100,'--','Color',[1 1 1]*0.5)
yticks([0 50 100])
grid on
ylabel('Probability of a test being positive','Interpreter','latex','FontSize',16)
ylim([-50,200])

legend([PlMsm Sh(1) Sh(4) PlPr],{
    'Measurements~~~'
    'Mean~~~'
    '95\% CI~~~'
    'Probability of a (wastewater) test being positive'
    },"Location","southeast",'Interpreter','latex','FontSize',16,'NumColumns',4)

exportgraphics(ax,"/home/ppolcz/Dropbox/Apps/Overleaf/COVID-22-Szennyviz-Revision1-Proba1/fig/GP_Cq_result.pdf","ContentType","vector")

%%
%}

if ~exist("Results/T_EMM.mat","file")
    GP_Exact_MM_Fast(hyp_Cp_Res)
    LogMean = Cp_mean*0;
    LogVar = Cp_var*0;
    LogCov = Cp_var*0;
    for i = 1:numel(Days)
        [LogMean(i),LogVar(i),LogCov(i)] = GP_Exact_MM_Fast(hyp_Cp_Res,Cp_mean(i),Cp_var(i));
        fprintf("%d/%d\n",i,numel(Days))
    end
    
    T_EMM = table(Days',Cp_mean,Cp_var,LogMean,LogVar,LogCov,'VariableNames',["Day","Cp_mean","Cp_var","LogMean","LogVar","LogCov"]);
    save("Results/T_EMM.mat","T_EMM")
else
    load("Results/T_EMM.mat")
end

T_Ms = join(T_Ms,T_EMM);

%{
%%
fig = figure(12);
fig.Position(3:4) = [1242 372];
delete(fig.Children);
ax = axes(fig);
hold on;
PlMsm = plot(T_Ms.Date,log10(T_Ms.PCR_Eredmeny),'o','MarkerSize',5,'Color',Color_1);
Sh = plot_mean_var(Start_Date+T_EMM.Day,T_EMM.LogMean,sqrt(T_EMM.LogVar),Color_2,"LineWidth",3);
grid on; box on

% Randoms = randn(numel(Days),1) .* sqrt(T_EMM.LogVar) + T_EMM.LogMean;
% PlResam = plot(Start_Date+Days,Randoms,'.','Color',Color_Black,'MarkerSize',10);
ax.XLim = Start_Date + [0 Days(end)];

Logger.latexify_axis(ax,14);

ylim([0 4])
ylabel('log10(PCR results)','Interpreter','latex','FontSize',14)

legend([ ...
    PlMsm Sh(1) Sh(4) ...
    ... PlResam ...
    ],{
    'Measurements (log10)'
    'Estimated mean'
    'Estimated 95\% CI'
    ... 'Resampled data'
    },"Location","northeast",'Interpreter','latex','FontSize',14)

title({'log10(PCR results): mean and variance computed from the GPs fitted to $k \mapsto C_\mathrm{q}$ and $C_\mathrm{q} \mapsto \mathrm{log10(PCR\ result)}$ using EMM'},'Interpreter','latex','FontSize',13)

exportgraphics(ax,"/home/ppolcz/Dropbox/Apps/Overleaf/COVID-22-Szennyviz-Revision1-Proba1/fig/GP_log10PCR_result.pdf","ContentType","vector")

%%
%}

% Log-normal distributions:
% https://en.wikipedia.org/wiki/Log-normal_distribution
% https://math.stackexchange.com/questions/660994/if-x-is-normal-is-expx-still-normal-how-to-find-its-mean-and-variance

%%

T_Ms.LoD = 1 ./ T_Ms.PCR_Bemert_RNS .* T_Ms.Vc ./ T_Ms.V0 .* T_Ms.Ve ./ T_Ms.Vs * 1e6;

Max_Nr = 150;
for i = 1:Max_Nr
    %% Ujrageneralni az egesz adatsort
    fprintf('%d/%d\n',i,Max_Nr)
    
    T = T_Ms;

    if i > 1
        n = height(T);
        
        % [1] Generate positivity rate --> rewrite .Pos
        T.Pos = rand(n,1) <= T.Pr_Pos;
        
        % [2] Generate Cp and LogResult
        e = cellfun(@(r,v1,c,v2) {r * chol([v1 c ; c v2])}, ...
            num2cell(randn(n,2),2), ...
            num2cell(T.Cp_var), ...
            num2cell(T.LogCov), ...
            num2cell(T.LogVar));
        e = vertcat(e{:});

        T.Cp = (T.Cp_mean + e(:,1)).* T.Pos;
        T.PCR_Eredmeny = 10.^(T.LogMean + e(:,2));
        T.Eredmeny_GcpL_szamitott_LOD = max(T.LoD,T.PCR_Eredmeny .* T.Pos .* T.Vc ./ T.V0 .* T.Ve ./ T.Vs * 1e6);
    end

    %% Remove outliers as proposed by 2021_Fernandez-Cassi.etal
    % Remove records with a higher/lower Cp than the median +- 1.5 cycles 
    
    T_Cp = table2timetable(unstack(T(:,["Date","Cp","Varos"]),"Cp","Varos","AggregationFunction",@max));
    % T_Cp_daily = retime(T_Cp,"daily");
    T_Cp_daily = sortrows(T_Cp);
    
    isout = isoutlier(T_Cp_daily,"movmedian",days(30)); % <----- RM OUTLIERS
    Vars = T_Cp_daily.Variables;
    Vars(isout == 1) = NaN;
    Vars = fillmissing(Vars,"nearest"); % <---- FILL MISSING VALUES
    Vars = movmean(Vars,10,1); % <------------- MOOVING MEDIAN FILTER (10 days)
    T_Cp_daily.Variables = Vars;
    
    T_medCp = stack(timetable2table(T_Cp_daily),T_Cp_daily.Properties.VariableNames, ...
        'IndexVariableName','Varos','NewDataVariableName','med_Cp');
    
    T = join(T,T_medCp,"Keys",["Date","Varos"]);
    %{
        % Ha a joint hibat dobna:
        Left = unique(T(:,["Date","Varos"]));
        Right = unique(T_medCp(:,["Date","Varos"]));
        ism = ismember(Left,Right);
        Left(~ism,:)
    %}
    
    % Find outliers, highly deviated measurements from the median Cq
    ldx = abs(T.Cp - T.med_Cp) > 3;
    
    %% Remove outliers
    T_nofilt = T;
    T_skipped = T(ldx,:);
    T(ldx,:) = [];
    
    %% Városonként lineárisan ábrázolni
    
    XLim = [
        datetime(2020,08,31)
        T.Date(end)
        ];
    
    FontSize = 12;
    Args = {"Interpreter","latex","FontSize",FontSize};
    
    % Csak ami kell egy külön táblába:
    R0 = T_nofilt(:,["Date","Varos","Eredmeny_GcpL_szamitott_LOD","Cp","LoD"]);
    R0 = renamevars(R0,"Eredmeny_GcpL_szamitott_LOD","Result");
    
    Varosok = unique(R0.Varos);
    
    % Csak ami kell egy külön táblába:
    R1 = T(:,["Date","Varos","Eredmeny_GcpL_szamitott_LOD","Cp","LoD"]);
    R1 = renamevars(R1,"Eredmeny_GcpL_szamitott_LOD","Result");
    
    % Csak ami kell egy külön táblába:
    Tmp = T_nofilt(:,["Date","Varos","Cp"]);
    T_Cp_daily_noft = table2timetable(unstack(Tmp,"Cp","Varos","AggregationFunction",@mean));
    
    % Limit of detection
    LoD = (1/2.5 * 30) * (1000 / 140) * (1000 / 50);
    
    Plot_Colors
    
    Col = 7;
    Row = 4;
    
    if i <= 0 % 2
        fig = figure(13120 + i);
        fig.Position(3:4) = [3583 1635];
        Tl = tiledlayout(Row*3,Col,'TileSpacing','compact','Padding','compact');
        % Tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
        
        for j = 1:min(numel(Varosok),prod(Tl.GridSize))
            varos = char(Varosok(j));
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
    
            Cp = T_Cp_daily_noft.(varos);
    
            ldx = ~isnan(Cp);
    
            plot(R0.Date(ldx0),R0.Cp(ldx0),'.-','MarkerSize',20,'Color',[1 1 1]*0.7)
            plot(R1.Date(ldx1),R1.Cp(ldx1),'.-','MarkerSize',20,'Color',Color_1)
            plot(T_Cp_daily.Date,T_Cp_daily.(varos),'LineWidth',3,'Color',Color_2)
            title("$C_\mathrm{q}$ " + VN_map(Varos),Args{:})
    
            grid on, box on
            xlim(XLim)
            ylim([0 40])
            if mod(i,Col) == 1
                ylabel('cycle (Cq)')
            end
    
            % drawnow
        end
    
        exportgraphics(fig,"/home/ppolcz/Dropbox/Apps/Overleaf/COVID-22-Szennyviz-Revision1-Proba1/fig/PCR-Results_All-" + num2str(i) + ".pdf","ContentType","vector")
    end
    
    %%
    
    XLim = [
        datetime(2020,08,31)
        T.Date(end)
        ];
    Date = XLim(1):XLim(2);
    
    FontSize = 12;
    Args = {"Interpreter","latex","FontSize",FontSize};
    

    if i <= 0 % 2
        
        for j = 1:numel(Varosok)
            varos = char(Varosok(j));
            %%
            varos = 'szeged';
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
            title("PCR results " + VN_map(Varos) + " \rule[1em]{0.001pt}{0.001pt}",Args{:});
            
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

            Cp = T_Cp_daily_noft.(varos);

            ldx = ~isnan(Cp);

            plot(R0.Date(ldx0),R0.Cp(ldx0),'.-','MarkerSize',20,'Color',[1 1 1]*0.7,'HandleVisibility','off')
            plot(R1.Date(ldx1),R1.Cp(ldx1),'.-','MarkerSize',20,'Color',Color_1,'HandleVisibility','off')
            plot(T_Cp_daily.Date,T_Cp_daily.(varos),'LineWidth',3,'Color',Color_2,'DisplayName','Moving median')
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

            tag = '';
            if isempty(ax.Title.String)
                tag = 'notit-';
            end
            exportgraphics(fig,"/home/ppolcz/Dropbox/Apps/Overleaf/COVID-22-Szennyviz-Revision1-Proba1/fig/PCR-Results_" + Varos + "-" + tag + "-" + num2str(i) + ".pdf","ContentType","vector")

            break
        end
    end    
    %%
    %%%
    
    % Csak ami kell egy külön táblába:
    Tmp = T(:,["Date","Eredmeny_GcpL_szamitott_LOD","Lakossag"]);
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
    for Tmp__ = 1:2
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
    
    if ~exist("Aggr","var")
        Aggr = Szv;
    else
        Aggr.Filt = [Aggr.Filt Szv.Filt];
    end
end

%%

% Plot_WData_City(R0,R1,T_Cp_daily_noft,'szeged');

%%

S = cov(Aggr.Filt(:,2:end)');

L = cholcov(S);
[n,m] = size(L);
Lp = [ L ; zeros(m-n,m) ];

Samples = Aggr.Filt(:,1) + Lp' * randn(size(Aggr.Filt(:,2:end)));

fig = figure(1231);
fig.Position(3:4) = [1112 255];
delete(fig.Children)
ax = axes(fig);
hold on, grid on, box on
% plot(Aggr.Date,Samples,'Color',[1 1 1]*0.6)
plot(Aggr.Date,round(Samples(:,4),-4),'LineWidth',3)
plot(Aggr.Date,Aggr.Filt(:,1),'LineWidth',3,'Color','red')

ylabel('Genome copy conc. [GC/L]','Interpreter','latex','FontSize',14)

xlim(Aggr.Date([1,end]))
ax.YLim(1) = 0;

Logger.latexify_axis(ax,12)

% exportgraphics(ax,"/home/ppolcz/Dropbox/Apps/Overleaf/COVID-22-Szennyviz-Revision1-Proba1/fig/Wastewater_curves.pdf","ContentType","vector")

%{
Aggr.Filt(:,2:end) = Samples;
save Results/Aggr.mat Aggr
%}