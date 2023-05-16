function [T,Quantity,Data] = Step5_ManyPRED_ODE(T,args)
arguments
    T
    args.Start_Date = datetime(2022,05,01);
    args.End_Date = T.Properties.UserData.Date_Last_Szennyviz + 90;
    args.ImLossRate = [0.005, 0.0055, 0.006, 0.0065, 0.007, 0.0075 ];
    args.RelTrRate = 1;
    args.TrRateMtp = (0.85:0.05:1.15).^1.5;
    % ---
    args.FigNr = 102;
    args.XLim = datetime(2022,01,01);
    args.TransitionLength = 30;
    args.TransitionLengthTr = [];
    args.TransitionLengthLs = [];
    args.TransitionLengthGn = [];
    args.TransitionLengthPr = [];
    % ---
    args.Plot = false;
end
%%
%  File: epid_tracking_NMPC.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study7_tracking
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2021. April 02. (2020b)
%  Minor revision on 2021. September 15. (2021a)
%  Revised on 2021. October 11. (2021b)
%  Revised on 2021. October 12. (2021b)
%  Major revision on 2022. January 07. (2021b)
% 

Date_Start_Pred = args.Start_Date;
Date_Last_H = T.Properties.UserData.Date_Last_Available_H;
Date_Last_Szv = T.Properties.UserData.Date_Last_Szennyviz;
Date_Last_Relevant_New_Cases = max(Date_Last_H - 4,Date_Last_Szv - 2);

R = T;

% Possible combinations
ls_TrRateMtp = args.TrRateMtp;
ls_PrH_Mtp = 0.9:0.1:1.1;

AllComb = array2table(allcomb(args.ImLossRate,ls_TrRateMtp,ls_PrH_Mtp),...
    "VariableNames",{'ImLossRate','TrRateMtp','PrH_Mtp'});

G_reset  

K = Epid_Par.GetK;

Date_End_PRED = args.End_Date;

Quantity = {

    % Predicted quantities
    @(T) T.ImLossRate , "Possible scenarios for the rate of immunity loss"
    @(T) T.TrRate , "Possible scenarios for the transmission rate"
    @(T) T.Param(:,K.pH), "Possible scenarios for the hospitalization probability"
    
    % Computed quantities
    @(T) T.New_Cases , "New cases"
    @(T) T.Infected , "All infected"
    @(T) T.H , "Hospitalized patients"
    
    };

Nr_Qty = height(Quantity);
Data = cell(Nr_Qty,height(AllComb));

%%

T_Prev = T(T.Date == Date_Start_Pred-1,:);
T = T(T.Date == Date_Start_Pred,:);

Np = T.Np;

D_TrRate = T.TrRate - T_Prev.TrRate;
D_ImLoss = T.ImLossRate - T_Prev.ImLossRate;
D_pH = T.Param(K.pH) - T_Prev.Param(K.pH);

T = tb_expand_timerange(T,[T.Date(1),Date_End_PRED]);
N = height(T);
for i = 1:N-1
    T.TrRate(i+1) = T.TrRate(i) + D_TrRate * (0.95)^i;
    T.ImLossRate(i+1) = T.ImLossRate(i) + D_ImLoss * (0.95)^i;
    T.Param(i+1,K.pH) = T.Param(i,K.pH) + D_pH * (0.95)^i;
end
T_clean = T(:,["S","L","P","I","A","H","D","R", "New_Cases", "Infected", "Infectious", "Rt", "Param", "TrRate", "ImLossRate", "ImGainRate", "TrRate_Ref", "L_iPeriod", "P_iPeriod", "I_iPeriod"]);

N_Tr = args.TransitionLengthTr; if isempty(N_Tr), N_Tr = args.TransitionLength; end
N_Ls = args.TransitionLengthLs; if isempty(N_Ls), N_Ls = args.TransitionLength; end
N_Gn = args.TransitionLengthGn; if isempty(N_Gn), N_Gn = args.TransitionLength; end
N_Pr = args.TransitionLengthPr; if isempty(N_Pr), N_Pr = args.TransitionLength; end

SigmoidTr = [Epid_Par.Sigmoid(1,0,N_Tr)' ; zeros(N-N_Tr,1)];
SigmoidLs = [Epid_Par.Sigmoid(1,0,N_Ls)' ; zeros(N-N_Ls,1)];
SigmoidGn = [Epid_Par.Sigmoid(1,0,N_Gn)' ; zeros(N-N_Gn,1)];
SigmoidPr = [Epid_Par.Sigmoid(1,0,N_Pr)' ; zeros(N-N_Pr,1)];

%%

[f,h,hvar,J] = epid_ode_model_8comp(Np);

for Scenario = 1:height(AllComb)
%%
    fprintf('Scenario %d / %d\n',Scenario,height(AllComb))
    T = T_clean;
    %%
        
    % Assumption 1. The stringency index is not bad.
    T.TrRate = SigmoidTr.*T.TrRate + (1-SigmoidTr).*(T.TrRate_Ref * args.RelTrRate * AllComb.TrRateMtp(Scenario));
    
    % Assumption 2. The immunity loss.
    T.ImLossRate = SigmoidLs.*T.ImLossRate + (1-SigmoidLs).*(T.ImLossRate*0 + AllComb.ImLossRate(Scenario));
    
    % Assumption 3. The immunity gain rate (by vaccination).
    ImGainRate_Future = 0;
    T.ImGainRate = SigmoidGn.*T.ImGainRate + (1-SigmoidGn).*ImGainRate_Future;
    
    % Assumption 4. Hospitalization probability.
    [Val,Idx] = min(T.Param(:,K.pH));
    pH_Multiplier = AllComb.PrH_Mtp(Scenario);
    T.Param(Idx:end,K.pH) = Val;
    T.Param(:,K.pH) = T.Param(:,K.pH) .* ((1 - pH_Multiplier)*SigmoidPr + pH_Multiplier);
    
    %{ 
    figure(623), stackedplot(T(:,["TrRate","ImLossRate","ImGainRate","Pr_H"]))
    %}
    
    %% Solve recursion -- the simplest way of prediction
        
    StateNames = string(cellfun(@char,num2cell(f.Input_1),'UniformOutput',false))';
    Tx = T(:,StateNames);
    
    x = Tx.Variables;
    for i = 1:N-1
        x(i+1,:) = x(i,:) + Fn_SLPIAHDR_ode(x(i,:)',T.Param(i,:)',T.TrRate(i),T.ImGainRate(i),T.ImLossRate(i))';
    end
    Tx.Variables = x;
    
    T(:,StateNames) = Tx;
    T = gen_depvars_from_LPIAHD(T);
    
    %% Collect data

    % H{Scenario} = T.H;
    % New_Cases{Scenario} = T.New_Cases;
    % Infected{Scenario} = T.Infected;

    for q = 1:Nr_Qty
        Data{q,Scenario} = Quantity{q,1}(T(1:end-1,:));
    end

end
%%

if ~args.Plot || args.FigNr <= 0
    return
end

load("../Data/T_Agens_2022_10_17.mat")

%%

if isscalar(args.XLim)
    args.XLim(2) = T.Date(end);
end

fig = figure(args.FigNr);
Tl = tiledlayout(3,2,'TileSpacing','tight','Padding','tight','TileIndexing','columnmajor');

ldx = R.Date <= Date_Last_Relevant_New_Cases;

Ax = [];
for q = 1:Nr_Qty
    Ax = [Ax ; nexttile];
    hold on; grid on; box on;

    qty = [Data{q,:}];
    Mean = mean(qty,2);
    Std = std(qty,0,2);

    if Quantity{q,2} == "New cases"
        plot(T_Agens.Date,T_Agens.New_Cases*5,'Color',[0.4660 0.6740 0.1880])
    elseif Quantity{q,2} == "Hospitalized patients"
        plot(T_Agens.Date,T_Agens.H,'Color',[0.4660 0.6740 0.1880])
    end        

    plot_mean_var(T.Date(1:end-1),Mean,Std,[0.9290 0.6940 0.1250],"FaceAlpha",0.2);
    % plot(T.Date(1:end-1),[Data{q,:}],'Color',[0.9290 0.6940 0.1250])
    plot(R.Date(ldx),Quantity{q,1}(R(ldx,:)),'LineWidth',2,'Color','red');
    ptitle(14,Quantity{q,2})
    drawnow
end

T.Properties.UserData.Link6 = linkprop(Ax,"XLim");
Ax(1).XLim = [R.Date(1) T.Date(end)];

DRange = R.Date(1):T.Date(end);
for ax = Ax'
    plot_vertical_milestone(ax,Date_Start_Pred,"Prediction");
    ax.XTick = DRange(day(DRange) == 1);
    ax.XMinorGrid = 'on';
    ax.XAxis.MinorTick = 'off';
    ax.XAxis.MinorTickValues = DRange(weekday(DRange) == 1);    
end

%%

fig = figure(args.FigNr+400);
Tl = tiledlayout(2,1,'TileSpacing','tight','Padding','tight','TileIndexing','columnmajor');

ldx = R.Date <= Date_Last_Relevant_New_Cases;

Ax = [];
for q = 1:Nr_Qty
    if ismember(Quantity{q,2},["New cases","Hospitalized patients"])
        Ax = [Ax ; nexttile];
        Logger.latexify_axis(Ax(end),12);
        hold on; grid on; box on;
    
        qty = [Data{q,:}];
        Mean = mean(qty,2);
        Std = std(qty,0,2);
    
        % if Quantity{q,2} == "New cases"
        %     PlAgent = plot(T_Agens.Date,T_Agens.New_Cases*5,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
        % elseif Quantity{q,2} == "Hospitalized patients"
        %     PlAgent = plot(T_Agens.Date,T_Agens.H,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
        % end        
    
        if Quantity{q,2} == "Hospitalized patients"
            PlOff = plot(R.Date,R.H_off_ma,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
        end

        Sh = plot_mean_var(T.Date(1:end-1),Mean,Std,[0.9290 0.6940 0.1250],"FaceAlpha",0.2);
        % plot(T.Date(1:end-1),[Data{q,:}],'Color',[0.9290 0.6940 0.1250])
        PlRec = plot(R.Date(ldx),Quantity{q,1}(R(ldx,:)),'LineWidth',2,'Color',[0 0.4470 0.7410]);
        title(Quantity{q,2},"Interpreter","latex","FontSize",16)
        drawnow
    end
end

T.Properties.UserData.Link6 = linkprop(Ax,"XLim");
Ax(1).XLim = args.XLim;

DRange = R.Date(1):T.Date(end);
for ax = Ax'
    plot_vertical_milestone(ax,Date_Start_Pred,"Prediction");
    ax.XTick = DRange(weekday(DRange) == 1);
    ax.XMinorGrid = 'on';
    ax.XAxis.MinorTick = 'off';
    ax.XAxis.MinorTickValues = DRange;
    
    legend(ax,[PlOff,PlRec,Sh(1),Sh(4)],["Official data","Reconstruction","Pred. -- mean","Pred. -- 95\% CI"],"Interpreter","latex","FontSize",12)
end


%%

fig = figure(args.FigNr+500);
Tl = tiledlayout(2,1,'TileSpacing','tight','Padding','tight','TileIndexing','columnmajor');

ldx = R.Date <= Date_Last_Relevant_New_Cases;

Ax = [];
for q = 1:Nr_Qty
    if ismember(Quantity{q,2},["Possible scenarios for the rate of immunity loss","Possible scenarios for the transmission rate"])
        Ax = [Ax ; nexttile];
        Logger.latexify_axis(Ax(end),12);
        hold on; grid on; box on;
    
        qty = [Data{q,:}];
        Mean = mean(qty,2);
        Std = std(qty,0,2);
        
        Sh = plot_mean_var(T.Date(1:end-1),Mean,Std,[0.9290 0.6940 0.1250],"FaceAlpha",0.2);
        PlRec = plot(R.Date(ldx),Quantity{q,1}(R(ldx,:)),'LineWidth',2,'Color',[0 0.4470 0.7410]);
        title(Quantity{q,2},"Interpreter","latex","FontSize",16)
        drawnow
    end
end

T.Properties.UserData.Link6 = linkprop(Ax,"XLim");
Ax(1).XLim = args.XLim;

DRange = R.Date(1):T.Date(end);
for ax = Ax'
    plot_vertical_milestone(ax,Date_Start_Pred,"Prediction");
    ax.XTick = DRange(weekday(DRange) == 1);
    ax.XMinorGrid = 'on';
    ax.XAxis.MinorTick = 'off';
    ax.XAxis.MinorTickValues = DRange;
    
    legend(ax,[PlRec,Sh(1),Sh(4)],["Reconstruction","Pred. -- mean","Pred. -- 95\% CI"],"Interpreter","latex","FontSize",12)
end

%%

fig = figure(args.FigNr+600);
Tl = tiledlayout(4,1,'TileSpacing','tight','Padding','tight','TileIndexing','columnmajor');

ldx = R.Date <= Date_Last_Relevant_New_Cases;

Ax = [];
for q = 1:Nr_Qty
    if ismember(Quantity{q,2},["New cases","Hospitalized patients","Possible scenarios for the rate of immunity loss","Possible scenarios for the transmission rate"])
        Ax = [Ax ; nexttile];
        Logger.latexify_axis(Ax(end),12);
        hold on; grid on; box on;
    
        qty = [Data{q,:}];
        Mean = mean(qty,2);
        Std = std(qty,0,2);
    
        % if Quantity{q,2} == "New cases"
        %     PlAgent = plot(T_Agens.Date,T_Agens.New_Cases*5,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
        % elseif Quantity{q,2} == "Hospitalized patients"
        %     PlAgent = plot(T_Agens.Date,T_Agens.H,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
        % end        
    
        if Quantity{q,2} == "Hospitalized patients"
            PlOff = plot(R.Date,R.H_off_ma,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
        end

        Sh = plot_mean_var(T.Date(1:end-1),Mean,Std,[0.9290 0.6940 0.1250],"FaceAlpha",0.2);
        % plot(T.Date(1:end-1),[Data{q,:}],'Color',[0.9290 0.6940 0.1250])
        PlRec = plot(R.Date(ldx),Quantity{q,1}(R(ldx,:)),'LineWidth',2,'Color',[0 0.4470 0.7410]);
        title(Quantity{q,2},"Interpreter","latex","FontSize",16)
        drawnow
    end
end
Ax(2).YLim(1) = 0;
Ax(1).YLim(1) = 0;

T.Properties.UserData.Link6 = linkprop(Ax,"XLim");
Ax(1).XLim = args.XLim;

DRange = R.Date(1):T.Date(end);
for ax = Ax'
    plot_vertical_milestone(ax,Date_Start_Pred,"Prediction");
    ax.XTick = DRange(weekday(DRange) == 1);
    ax.XMinorGrid = 'on';
    ax.XAxis.MinorTick = 'off';
    ax.XAxis.MinorTickValues = DRange;

    if ax ~= Ax(end)
        ax.XTickLabel = [];
    end
end

legend(Ax(end),[PlOff,PlRec,Sh(1),Sh(4)],["Official data","Reconstruction","Pred. -- mean","Pred. -- 95\% CI"],"Interpreter","latex","FontSize",12,"Location","northwest")

TODAY = datestr(date,29);
DIR = "/home/ppolcz/T/Dropbox/Peti/Munka/01_PPKE_2020/COVID-Elorejelzesek/UIO_and_Opt/UIO_and_Opt_" + TODAY;
if ~exist(DIR,"dir")
    mkdir(DIR);
end
exportgraphics(fig,DIR + "/Pred_" + TODAY + ".pdf");

Ax(1).XLim(1) = datetime(2021,10,24);
exportgraphics(fig,DIR + "/Pred_" + TODAY + "_From_2021Okt24.pdf");

Ax(1).XLim(1) = datetime(2020,08,20);
Ax(4).YLim(2) = 13000;
exportgraphics(fig,DIR + "/Pred_" + TODAY + "_From_2020Aug20.pdf");

end
