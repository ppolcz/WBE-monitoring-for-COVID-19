%%
%  File: Main.m
%  Directory: /home/ppolcz/T/_Epid/Approach_UIO_then_Opt/Ver_2022_05_31_UIO_then_Opt
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. June 07. (2022a)
%
% References
% 
%  [1] B. Csutak, P. Polcz and G. Szederkényi, "Model-based epidemic data
%  reconstruction using feedback linearization," 2022 International
%  Conference on Electrical, Computer and Energy Technologies (ICECET),
%  2022, pp. 1-6, doi: 10.1109/ICECET55527.2022.9873061.

G_reset
Plot_Colors

%% Load parameters from xls

% Read the appropriate parameter table (.xlsx), then, generate the
% parameter trajectories
[P,Q,Params_xls] = Epid_Par.Get("HUN_Reconstruction");
Q_WW = Epid_Par.Read("HUN_Wastewater_Mtp");

% Helper table
Q_nc = Q(:,["Date","TransitionLength","New_Cases_Shift","New_Cases_Mtp"]);
Q_nc = renamevars(Q_nc,["New_Cases_Mtp","New_Cases_Shift"],["Mtp","Shift"]);

% Load all data
[T_Data,~,s] = Step0_Preparing_Data(P);
UserData = T_Data.Properties.UserData;

P(P.Date < T_Data.Date(1) | T_Data.Date(end) < P.Date,:) = [];
T_Data.Param = P.Param;

%% Compute initial estimates for L,P,I,A,H,D using an unknown input observer

DIR = "Results";
if ~exist(DIR,"dir")
    mkdir(DIR)
end
xls = DIR + "/T_1_UIO.xls";
if ~exist(xls,"file")
    % Dynamic inversion using hospitalization data only
    T_UIO = Step1_LTI_inversion_UIO(T_Data,s);
    
    writetimetable(T_UIO,xls,"Sheet","T");
    writetable(Q,xls,"Sheet","Q","WriteRowNames",true);
else
    T_UIO = pcz_read_epid_table(UserData,xls,"Sheet","T");
end

%% Calibrate L,P,I,A,H,D through optimization using wastewater data

xls = DIR + "/T_2_Opt.xls";
if ~exist(xls,"file")

    % Reconstruction and parameter calibration using wastewater
    [T,Q_Opt] = Step2_LTI_inversion_Opt(T_UIO,Q,Q_WW,"UseWasteWater",true);
    
    % Detect multipliers for detected cases
    [~,~,~,~,Q_nc] = pcz_detect_multiplier(Q_nc,T.Date,T.New_Cases,T.New_Cases_off,"Modify",true);
    Q_Opt.New_Cases_Mtp = Q_nc.Mtp;

    writetimetable(T,xls,"Sheet","T");
    writetable(Q_Opt,xls,"Sheet","Q","WriteRowNames",true);
else
    T = pcz_read_epid_table(UserData,xls,"Sheet","T");
    Q_Opt = readtable(xls,"Sheet","Q","ReadRowNames",true);
end

%% Full reconstruction using an iterative approach

% Enhanced vaccination model to compute waning and susceptibles
T_Full = Step3_Opt_Iterative(T);

xls = DIR + "/T_3_Full.xls";
writetimetable(T_Full,xls,"Sheet","T");
writetable(Q_Opt,xls,"Sheet","Q","WriteRowNames",true);


%% Compute variances when parameters are uncertain

T_Full = Step4_Parameter_uncertainty(T_Full);
Plot_Reconstruction_MeanStd(T_Full,Q)

xls = DIR + "/T_4_Full_ParUnc.xls";
writetimetable(T_Full,xls,"Sheet","T");
writetable(Q_Opt,xls,"Sheet","Q","WriteRowNames",true);

%% Visualize

% 3.3 Detect multipliers for detected cases
Q_nc = Q(:,["Date","TransitionLength","New_Cases_Shift","New_Cases_Mtp"]);
Q_nc = renamevars(Q_nc,["New_Cases_Mtp","New_Cases_Shift"],["Mtp","Shift"]);
[New_Cases_Off_Shifted,LegEntries_NC,~,~,Q_nc] = pcz_detect_multiplier(Q_nc,T_Full.Date,T_Full.New_Cases,T_Full.New_Cases_off,"Modify",true);
Q.New_Cases_Mtp = Q_nc.Mtp;

Q_nc_Sh0 = Q_nc;
Q_nc_Sh0.Shift = Q_nc.Shift*0;
[New_Cases_Off,LegEntries_NC,~] = pcz_detect_multiplier(Q_nc_Sh0,T_Full.Date,T_Full.New_Cases,T_Full.New_Cases_off, ...
    "Modify",false,"RemoveLast",1,"Legend","Detected cases");

Plot_Reconstruction(T_Full,T_Full,Q_Opt,New_Cases_Off,LegEntries_NC,Q_nc);

%% Compute statistics

Qil = Q_WW;
Qil('Alpha',:).Date = datetime(2021,01,01);
Qil('Delta',:).Date = datetime(2021,04,24);
Qil('OmicronBA12',:).Date = datetime(2021,12,15);
Qil('OmicronBA5',:).Date = datetime(2022,05,01);
Qil('OmicronBQ1',:).Date = datetime(2022,11,01);

Np = T_Full.Np(1);

T = T_Full(1:end-1,:);
P = epid_get_variant_dominance_pattern(Q_WW,T.Date([1,end]));
Pil = epid_get_variant_dominance_pattern(Qil,T.Date([1,end]));

Pattern = [
    P.V_Original ...
    P.V_Alpha ...
    P.V_Delta ...
    P.V_OmicronBA12 ...
    P.V_OmicronBA5 ...
    P.V_OmicronBQ1
    ];

PatternIl = [
    Pil.V_Original ...
    Pil.V_Alpha ...
    Pil.V_Delta ...
    Pil.V_OmicronBA12 ...
    Pil.V_OmicronBA5 ...
    Pil.V_OmicronBQ1 ...
    ];

S = table('RowNames',[Qil.Properties.RowNames(2:end-1) ; 'SUM']);

S.Wave_Start = [ Q_WW.Date(2:end-1) ; Q_WW.Date(2) ];
S.Wave_End = [ Q_WW.Date(3:end) ; Q_WW.Date(end) ];

for i = 1:height(S)-1
    S.All_Cases(i) = round(sum(T.New_Cases .* Pattern(:,i)));
    S.Perc_Cases(i) = round(100 * S.All_Cases(i) / Np);
end

S.ImLoss_Start = [ Qil.Date(2:end-1) ; Qil.Date(2) ];
S.ImLoss_End = [ Qil.Date(3:end) ; Qil.Date(end) ];

for i = 1:height(S)-1
    S.All_ILoss(i) = round(sum(T.New_ImLoss .* PatternIl(:,i)));
    S.All_IGain(i) = round(sum(T.New_ImGain .* PatternIl(:,i)));
    S.All_Vacc(i) = round(sum(diff(T_Full.VacF + T_Full.VacB) .* PatternIl(:,i)));

    S.Perc_ILoss(i) = round(100 * S.All_ILoss(i) / Np);
    S.Perc_IGain(i) = round(100 * S.All_IGain(i) / Np);
    S.Perc_Vacc(i) = round(100 * S.All_Vacc(i) / Np);

    S.Eff_Vacc(i) = round(100 * S.All_IGain(i) / S.All_Vacc(i));

    S.Daily_ILoss(i) = round(S.All_ILoss(i) / sum(PatternIl(:,i)));
end

% Numeric columns of table S
Idx = find(varfun(@isnumeric,S,'output','uniform'));
for i = Idx
    S(end,i) = table(sum(S(1:end-1,i).Variables,1,'omitnan'));
end

xls = DIR + "/Results.xls";
writetable(S,xls,"Sheet","Q","WriteRowNames",true);

%% Long-term prediction

Date_Start_Pred = datetime(2022,05,01);
[T_PRED,Qty,Dat2] = Step5_ManyPRED_ODE(T_Full,"Start_Date",Date_Start_Pred,"FigNr",4014,...
    "TransitionLengthLs",90);

Plot_Prediction(T_Full,T_PRED,Q_Opt,Qty,Dat2, ...
    "Date_Start_Pred", Date_Start_Pred);

