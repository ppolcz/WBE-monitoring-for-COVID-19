function [ret] = Plot_Reconstruction_MeanStd(T,Q,args)
arguments
    T,Q
    args.FontSize = 13;
    args.FigNr = 2311;
end
%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2022. October 21. (2022b)
%

Plot_Colors

%% Visualization

XLim = [ datetime(2020,08,01) T.Properties.UserData.Date_Last_Szennyviz ];

fig = figure(args.FigNr);
% fig.Position(3) = 1494;
% fig.Position(4) = 1000;

Tl = tiledlayout(15,1,"TileSpacing","tight","Padding","compact");

% ax = subplot(511); hold on, plot_mean_var(T.Date,T.ImLossRate,T.ImLossRateStd), 
% ax = subplot(512); hold on, plot_mean_var(T.Date,T.TrRate,T.TrRateStd), ylim([0,3])
% ax = subplot(513); hold on, plot_mean_var(T.Date,T.L,T.StateStd(:,2)), ax.YLim(1) = 0;
% ax = subplot(514); hold on, plot_mean_var(T.Date,T.New_Cases,T.OutputStd(:,2)), ax.YLim(1) = 0;
% ax = subplot(515); hold on, plot_mean_var(T.Date,T.H,real(T.StateStd(:,J.H))), ax.YLim(1) = 0;

tidx = 1;
% -------
[Ax,tidx] = nt(tidx,5);
plot_mean_var(T.Date,T.TrRate,T.TrRateStd);
ylim([0,3])

ptitle(12,'Transmission rate')

% Leg = legend([Bar Pl Pl2],strsplit(LegEntries_NC + '; New cases (with wastewater); New cases (without wastewater)',"; "), ...
%     "Interpreter","latex","FontSize",args.FontSize,"Location","northwest","Box","on");
% Leg.NumColumns = 2;

[Ax(end+1),tidx] = nt(tidx,5);
plot_mean_var(T.Date,T.ImLossRate,T.ImLossRateStd);
ylim([0,0.02])

ptitle(12,'Immunity loss rate')

% Leg(end+1) = legend( ...
%     'Rate of immunity loss ($\omega$)', ...
%     'Immunization rate ($\nu$)', ...
%     "Interpreter","latex","FontSize",args.FontSize,"Location","northwest");
% ptitle(12,'Estimated waning (proportion of \textbf{R} loosing immunity within 24 hours) and immunization (proportion of \textbf{S} gaining immunity within 24 hours)')

[Ax(end+1),tidx] = nt(tidx,5);
plot_mean_var(T.Date,T.New_Cases,T.OutputStd(:,2));
Ax(end).YLim(1) = 0;

ptitle(12,'New cases')

% -------

for idx = 1:numel(Ax)
    ax = Ax(idx);
    Logger.latexify_axis(ax,args.FontSize);
    % ax.XTick = T.Date(day(T.Date) == 1 & mod(month(T.Date),3) == 1);
    ax.XTick = T.Date(day(T.Date) == 1 & mod(month(T.Date),3) == 0);
    ax.XMinorGrid = 'on';
    ax.XAxis.MinorTick = 'off';
    ax.XAxis.MinorTickValues = T.Date(day(T.Date) == 1);    

    ax.Parent.UserData.CompName = "Pred_MPC";

    ax.XLim = XLim;
    
    % if ismember(idx,[1,7])
    %     plot_vertical_variants(ax,Q,"Except",["School","TurnPoint2"])
    % else
        Plot_Vertical_variants(ax,Q,"Except",["School","TurnPoint2","[VER]New variant"])
    % end
    drawnow 
end

for i = 1:numel(Leg)
    Leg(i).Position(2) = Leg(i).Position(2) - 0.005;
end

return
%%
fig = gcf;
fig.Position(3) = 1494;
fig.Position(4) = 1825;

Today = datetime("today");
Today.Format = "uuuu-MM-dd";
fMain = "/home/ppolcz/Dropbox/Peti/Munka/01_PPKE_2020/Dokumentaciok/Docs_CsutakB_PhD/07_COVID-Szennyviz/actual/fig/Reconstruction-" + string(Today) + ".pdf";
clipboard("copy",fMain)
exportgraphics(fig,fMain,'ContentType','vector')

end

function [ax,idx] = nt(idx,w)
    ax = nexttile(idx,[w 1]);
    hold on, grid on, box on
    idx = idx + w;
end