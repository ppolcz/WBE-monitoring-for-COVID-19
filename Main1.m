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

NrValsOneSide = 1;
NrVals = 2*NrValsOneSide + 1;

ls = logspace(0,2,NrValsOneSide);
ls = [-ls(end:-1:1) ls] / 100;

stdperc = {
    "Period_L"    10
    "Period_P"    10
    "Period_A"    10
    "Period_I"    10
    "Period_H"    10
    "Rel_beta_A"  10
    "Pr_I"        10
    "Pr_H"        10
%     "Pr_D"        10
    }.';
stdperc = struct(stdperc{:});

xls_fh = @(fn,j) DIR + "/T__" + fn + "__" + sprintf("%02d",j) + ".xls";
fig_fh = @(fn,qty) DIR + "/Fig__" + qty + "__vs__" + fn + ".png";

mtpy = struct;
for fn = string(fieldnames(stdperc))'
    Qi = Q;

    mtpy.(fn) = 1 + ls * stdperc.(fn) / 100;
    vals = Q.(fn) * mtpy.(fn);

%     for j = 1:NrVals-1
%         
%         Qi.(fn) = vals(:,j);
%         Pi = Epid_Par.Generate_Timetable(Qi,Params_xls);
% 
%         Ti = T_Data;
%         Ti.Param = Pi.Param(Ti.Date(1) <= Pi.Date & Pi.Date <= Ti.Date(end),:);
%     
%         Ti = Step1_LTI_inversion_UIO(Ti,s);
%         [Ti,Qi] = Step2_LTI_inversion_Opt(Ti,Qi,"UseWasteWater",true);
%     
%         [~,~,~,~,Q_nc] = pcz_detect_multiplier(Q_nc,Ti.Date,Ti.New_Cases,Ti.New_Cases_off,"Modify",true);
%         Qi.New_Cases_Mtp = Q_nc.Mtp;
%         Qi.qC = 1 ./ Qi.Szennyviz_Mtp;
%         
%         xls = xls_fh(fn,j);
%         writetimetable(Ti,xls,"Sheet","T");
%         writetable(Qi,xls,"Sheet","Q");
%     end
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
    Qty_T.(qty) = zeros(height(T),NrVals);
    Qty_T.(qty)(:,NrValsOneSide+1) = T.(qty);
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
    Qty_Q.(qty) = zeros(height(Q_Opt),NrVals);
    Qty_Q.(qty)(:,NrValsOneSide+1) = Q_Opt.(qty);
end

ColNames = [Quantities_T,Quantities_Q];
RowNames = fieldnames(stdperc);
Variables = num2cell(nan(numel(RowNames),numel(ColNames)),1);
Tbl = table(Variables{:},'VariableNames',ColNames,'RowNames',RowNames);

Idx = [ 1:(NrVals-1)/2 , (NrVals+1)/2+1:NrVals ];
for fn = string(fieldnames(stdperc))'
    m = ones(1,NrVals);
    m(Idx) = mtpy.(fn);
    LegEntries = cellfun(@(a) sprintf("%s: %5.1f%%",strrep(fn,"_"," "),a*100),num2cell(m));
        
    for j = 1:NrVals-1
        xls = xls_fh(fn,j);
        fprintf('Loading %s ...\n',xls);
        Ti = pcz_read_epid_table(UserData,xls,"Sheet","T");
        Qi = readtable(xls,"Sheet","Q");

        [~,~,~,~,Q_nc] = pcz_detect_multiplier(Q_nc,Ti.Date,Ti.New_Cases,Ti.New_Cases_off,"Modify",true);
        Qi.New_Cases_Mtp = Q_nc.Mtp;
        writetable(Qi,xls,"Sheet","Q")

        % Collect variables from T
        for qty = Quantities_T
            Qty_T.(qty)(:,Idx(j)) = Ti.(qty);
        end

        % Collect variables from Q
        for qty = Quantities_Q
            Qty_Q.(qty)(:,Idx(j)) = Qi.(qty);
        end
    end

    Qty_Q

    Norm = @(x) sum(abs(x));
    for qty = Quantities_T
        dQty = Qty_T.(qty)(1:end-1,:) - T.(qty)(1:end-1,:);
        mQty = T.(qty)(1:end-1,:);

        norm_Qty = max(Norm(dQty) / Norm(mQty));
        Tbl(fn,qty) = table(norm_Qty);
    end
    for qty = Quantities_Q
        dQty = Qty_Q.(qty) - Q_Opt.(qty);
        mQty = Q_Opt.(qty);

        norm_Qty = max( Norm(dQty) / Norm(mQty) );
        Tbl(fn,qty) = table(norm_Qty);
    end

    continue

    fig = figure(1231);
    fig.Position(3:4) = [1920,500];

    for qty = Quantities_T
        plot(T.Date,Qty_T.(qty)), legend(LegEntries,"Location","northwest");
        title("Epidemic curve [" + strrep(qty,"_"," ") + "] reconstructed for different [" + strrep(fn,"_"," ") + "]")
        ax = gca;
        ax.YLim(1) = 0;
        ax.XLim = T.Date([1,end]);
        grid on
        switch qty
            case "Rc"
                ax.YLim(2) = 3;
        end
        drawnow
    
        exportgraphics(fig,fig_fh(fn,qty))
    end

    for qty = Quantities_Q
        plot(Qty_Q.(qty),'.-','MarkerSize',20), legend(LegEntries,"Location","northwest");
        title("Parameter [" + strrep(qty,"_"," ") + "] estimated for different [" + strrep(fn,"_"," ") + "]")
        ax = gca;
        ax.XLim = [1-0.5,height(Qty_Q)+0.5];
        ax.XTick = 1:height(Qty_Q);
        ax.XTickLabel = Qty_Q.Label;
        grid on
        drawnow
    
        exportgraphics(fig,fig_fh(fn,qty))
    end
end

Tbl.Variables = round(Tbl.Variables * 100);
Tbl(["Pr_H"],:) = [];

%%

% set(groot,'defaulttextinterpreter','latex');  
% set(groot, 'defaultAxesTickLabelInterpreter','latex');  
% set(groot, 'defaultLegendInterpreter','latex'); 

RowNames_LaTeX = {
    '$\tau_{\mathrm{L}}$'
    '$\tau_{\mathrm{P}}$'
    '$\tau_{\mathrm{A}}$'
    '$\tau_{\mathrm{I}}$'
    '$\tau_{\mathrm{H}}$'
    '$q_{\mathrm{A}}$'
    '$p_{\mathrm{I}}$'
    % '$p_{\mathrm{H}}$'
    };

ColNames_LaTeX = {
    '$u$'
    '$\mathbf{K}\!-\!\mathbf{D}$'
    '$R_{\mathrm{c}}$'
    '$p_{\mathrm{H}}$'
    '$p_{\mathrm{D}}$'
    '$q_{\mathrm{c}}$'
    '$q_{\mathrm{u}}$'
    };

fig = figure(45);
fig.Position(3:4) = [557 305];
H = heatmap(ColNames_LaTeX,RowNames_LaTeX,Tbl.Variables);
H.XLabel = "Computed quantities";
H.YLabel = "Modified parameters";

H.NodeChildren(3).FontSize = 12;
H.NodeChildren(3).TickLabelInterpreter = "latex";
H.NodeChildren(3).XAxis.Label.Interpreter = "latex";
H.NodeChildren(3).YAxis.Label.Interpreter = "latex";

% Weight values in colorbar
H.NodeChildren(2).FontSize = 12;
H.NodeChildren(2).TickLabelInterpreter = "latex";
H.NodeChildren(2).Label.Position(1) = 3;
H.NodeChildren(2).Label.Interpreter = "latex";
H.NodeChildren(2).Label.String = "Relative distance in \%";
H.NodeChildren(2).Label.Rotation = 270;
H.NodeChildren(2).Label.FontSize = 12;

H.NodeChildren(1).TickLabelInterpreter = "latex";

exportgraphics(fig,"/home/ppolcz/Dropbox/Apps/Overleaf/COVID-22-Szennyviz-Revision1-Proba1/fig/Parameter_sensitivity.pdf","ContentType","vector")
