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

% Helper table
Q_nc = Q(:,["Date","TransitionLength","New_Cases_Shift","New_Cases_Mtp"]);
Q_nc = renamevars(Q_nc,["New_Cases_Mtp","New_Cases_Shift"],["Mtp","Shift"]);

% Load all data
[T_Data,~,s,~] = Step0_Preparing_Data(P);
UserData = T_Data.Properties.UserData;

P(P.Date < T_Data.Date(1) | T_Data.Date(end) < P.Date,:) = [];
T_Data.Param = P.Param;

DIR = "/home/ppolcz/Dropbox/Peti/Munka/01_PPKE_2020/COVID-Elorejelzesek/UIO_and_Opt/Sensitivity_2023-04-24";
DIR = "Results";
xls = DIR + "/T_UIO.xls";
if ~exist(xls,"file")
    % Dynamic inversion using hospitalization data only
    T_UIO = Step1_LTI_inversion_UIO(T_Data,s);
    
    writetimetable(T_UIO,xls,"Sheet","T");
    writetable(Q,xls,"Sheet","Q");
else
    T_UIO = pcz_read_epid_table(UserData,xls,"Sheet","T");
end

xls = DIR + "/T.xls";
if ~exist(xls,"file")

    % Reconstruction and parameter calibration using wastewater
    [T,Q_Opt] = Step2_LTI_inversion_Opt(T_UIO,Q,"UseWasteWater",true);
    
    % Detect multipliers for detected cases
    [~,~,~,~,Q_nc] = pcz_detect_multiplier(Q_nc,T.Date,T.New_Cases,T.New_Cases_off,"Modify",true);
    Q_Opt.New_Cases_Mtp = Q_nc.Mtp;
    Q_Opt.qC = 1 ./ Q_Opt.Szennyviz_Mtp;

    writetimetable(T,xls,"Sheet","T");
    writetable(Q_Opt,xls,"Sheet","Q");
else
    T = pcz_read_epid_table(UserData,xls,"Sheet","T");
    Q_Opt = readtable(xls,"Sheet","Q");
end


if false
    %% Full reconstruction
    
    % Enhanced vaccination model to compute waning and susceptibles
    T_Full = Step3_Opt_Iterative(T);
    
    T_Full = Step4_Parameter_uncertainty(T_Full);
    Plot_Reconstruction_MeanStd(T_Full,Q)
    
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
end

%%

load Results/Aggr.mat

NrSim = width(Aggr.Filt);

xls_fh = @(j) DIR + "/T__Szv__" + sprintf("%03d",j) + ".xls";
fig_fh = @(qty) DIR + "/Fig__" + qty + "__vs__Szv.png";

return

for i = 1:NrSim-1
    Ti = T_UIO;
    Idx_Szv = Aggr.Date(1) <= Ti.Date & Ti.Date <= Aggr.Date(end);
    Ti.Szennyviz(Idx_Szv) = max(Aggr.Filt(:,i+1),0);

    [Ti,Qi] = Step2_LTI_inversion_Opt(Ti,Q,"UseWasteWater",true);

    [~,~,~,~,Q_nc] = pcz_detect_multiplier(Q_nc,Ti.Date,Ti.New_Cases,Ti.New_Cases_off,"Modify",true);
    Qi.New_Cases_Mtp = Q_nc.Mtp;
    Qi.qC = 1 ./ Qi.Szennyviz_Mtp;
    
    xls = xls_fh(i);
    writetimetable(Ti,xls,"Sheet","T");
    writetable(Qi,xls,"Sheet","Q");
end

%%

% Collected variables from T
Quantities_T = [
    "New_Cases"
    "Infected"
    "Rc"
    ].';
Date = T.Date;
Qty_T = timetable(Date);
for qty = Quantities_T
    Qty_T.(qty) = nan(height(T),NrSim);
    Qty_T.(qty)(:,1) = T.(qty);
end

% Collected variables from Q
Quantities_Q = [
    "Pr_H"
    "Pr_D"
    "Szennyviz_Mtp"
    "New_Cases_Mtp"
    ].';
Qty_Q = Q(:,["Date","TransitionLength"]);
for qty = Quantities_Q
    Qty_Q.(qty) = nan(height(Q_Opt),NrSim);
    Qty_Q.(qty)(:,1) = Q_Opt.(qty);
end

for i = 1:NrSim-1
    xls = xls_fh(i);
    if ~exist(xls,"file")
        Aggr.Filt(:,i+1) = Aggr.Filt(:,i+1) * NaN;
        continue
    end

    fprintf('Loading %s ...\n',xls);
    Ti = pcz_read_epid_table(UserData,xls,"Sheet","T");
    Qi = readtable(xls,"Sheet","Q");

    % Collect variables from T
    for qty = Quantities_T
        Qty_T.(qty)(:,i+1) = Ti.(qty);
    end

    % Collect variables from Q
    for qty = Quantities_Q
        Qty_Q.(qty)(:,i+1) = Qi.(qty);
    end
end

%%

fig = figure(1231);
fig.Position(3:4) = [1112 355];
delete(fig.Children)

Tl = tiledlayout(1,1);
ax = nexttile;
hold on, grid on, box on
plot(Aggr.Date,Aggr.Filt(:,2:end),'Color',[1 1 1]*0.6)
plot(Aggr.Date,Aggr.Filt(:,1),'LineWidth',3,'Color','red')
% --
ylabel('Genome copy conc. [GC/L]','Interpreter','latex','FontSize',14)
Plot_Vertical_variants(ax,Q,"Except",["Wild"])
Logger.latexify_axis(ax,12)
% --
xlim(Aggr.Date([1,end]))
ax.YLim(1) = 0;

ax = nexttile("south");
hold on, grid on, box on
plot(T.Date,T.H_off_ma,'LineWidth',2);
% --
ylabel('$\mathbf{H}^{\mathrm{Off}}$','Interpreter','latex','FontSize',14)
Plot_Vertical_variants(ax,Q,"Except",["Wild"])
Logger.latexify_axis(ax,12)
% -- 
xlim(Aggr.Date([1,end]))

exportgraphics(fig,"/home/ppolcz/Dropbox/Apps/Overleaf/COVID-22-Szennyviz-Revision1-Proba1/fig/Wastewater_curves.pdf","ContentType","vector")


%%

Qty_MeanT = Qty_T(:,[]);
Qty_StdT = Qty_T(:,[]);

Qty_MeanQ = Qty_Q(:,[]);
Qty_StdQ = Qty_Q(:,[]);

for qty = Quantities_T
    Qty_MeanT.(qty) = mean(Qty_T.(qty),2,'omitnan');
    Qty_StdT.(qty) = std(Qty_T.(qty),0,2,'omitnan');
end

for qty = Quantities_Q
    Qty_MeanQ.(qty) = mean(Qty_Q.(qty),2,'omitnan');
    Qty_StdQ.(qty) = std(Qty_Q.(qty),0,2,'omitnan');
end

Qty_CIpercQ = Qty_StdQ;
Qty_CIpercQ.Variables = round(Qty_StdQ.Variables * 2 ./ Qty_MeanQ.Variables * 100);

Qty_Q(["Transient","TurnPoint2","Future"],:) = [];
Qty_MeanQ(["Transient","TurnPoint2","Future"],:) = [];
Qty_CIpercQ(["Transient","TurnPoint2","Future"],:) = [];
writetable(Qty_MeanQ,"Results/Qty_CIpercQ.xls","Sheet","Mean")
writetable(Qty_CIpercQ,"Results/Qty_CIpercQ.xls","Sheet","CI")

%%

fig = figure(1231);
fig.Position(3:4) = [1920,500];

LegEntries = cellfun(@(i) {num2str(i)},num2cell(0:NrSim-1));

for qty = Quantities_T
%%
    plot(T.Date,Qty_T.(qty)); 
    % legend(LegEntries{:});
    title("Epidemic curve [" + strrep(qty,"_"," ") + "] reconstructed for different wastewater curves")
    ax = gca;
    ax.YLim(1) = 0;
    ax.XLim = T.Date([1,end]);
    grid on
    switch qty
        case "Rc"
            ax.YLim(2) = 3;
    end
    drawnow

    exportgraphics(fig,fig_fh(qty))
end

for qty = Quantities_Q
    plot(Qty_Q.(qty),'.-','MarkerSize',20), 
    % legend(LegEntries,"Location","northwest");
    title("Parameter [" + strrep(qty,"_"," ") + "] estimated for different wastewater curves")
    ax = gca;
    ax.XLim = [1-0.5,height(Qty_Q)+0.5];
    ax.XTick = 1:height(Qty_Q);
    ax.XTickLabel = Qty_Q.Label;
    grid on
    drawnow

    exportgraphics(fig,fig_fh(qty))
end


hold off
plot(Aggr.Date,Aggr.Filt(:,2:end),'Color',[1 1 1]*0.6)
hold on
plot(Aggr.Date,Aggr.Filt(:,1),'LineWidth',3,'Color','red')
grid on

ax = gca;
ax.XLim = Aggr.Date([1,end])
ax.YLim(1) = 0;
grid on
drawnow

exportgraphics(fig,"Results/Fig__Aggr.png")


%%

% set(groot,'defaulttextinterpreter','latex');  
% set(groot, 'defaultAxesTickLabelInterpreter','latex');  
% set(groot, 'defaultLegendInterpreter','latex'); 

