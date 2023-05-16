function [ret] = Plot_Reconstruction(T0,T,Q,New_Cases_Off,LegEntries_NC,Q_nc,args)
arguments
    T0,T,Q,New_Cases_Off,LegEntries_NC,Q_nc
    args.FontSize = 13;
    args.FigNr = 1311;
end
%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2022. October 21. (2022b)
%

Plot_Colors

%% Visualization

if exist("T_Full_H","var") && exist("T_Full_SzH","var") && exist("Q_Opt","var") ...
    && exist("New_Cases_Off","var") && exist("LegEntries_NC","var") && exist("Q_nc","var") ...
    && exist("T_Ms","var") 
    Plot_Reconstruction(T_Full_H,T_Full_SzH,Q_Opt,New_Cases_Off,LegEntries_NC,Q_nc);
    return
end

[~,~,~,J] = epid_ode_model_8comp;

XLim = [ datetime(2020,08,01) T.Properties.UserData.Date_Last_Szennyviz ];

fig = figure(args.FigNr);
fig.Position(3) = 1494;
fig.Position(4) = 1825;

Tl = tiledlayout(30,1,"TileSpacing","tight","Padding","compact");

tidx = 1;
% -------
[Ax,tidx] = nt(tidx,5);
Bar = bar(T.Date,New_Cases_Off',1,'stacked','FaceAlpha',0.7);
for j = 1:numel(Bar)
    Bar(j).DisplayName = sprintf("~~\\makebox[1cm]{$\\times\\, %g$\\hfill} (" + New_Cases_Colors{j,2} + ")",Q_nc.Mtp(j));
    Bar(j).FaceColor = New_Cases_Colors{j,1};    
end

% Pl = plot(T.Date,T.New_Cases,'k');
Sh = plot_mean_var(T.Date,T.New_Cases,T.OutputStd(:,J.Daily_New),[0 0 0],"PMLineStyle",'-',"LineWidth",1);
Pl = Sh(1);
Pl2 = plot(T0.Date,T0.New_Cases,'k--');
Pl_Unc = area(datetime(1991,03,20) + [0 1],[1 1],'FaceColor',[0 0 0],'FaceAlpha',0.2,'EdgeColor',[0,0,0]);
Leg = legend([Bar Pl Pl_Unc Pl2],strsplit(LegEntries_NC + '; New cases with wastewater; and its uncertainty ($\pm 95\%$ CI); New cases (without wastewater)',"; "), ...
    "Interpreter","latex","FontSize",args.FontSize,"Location","northwest","Box","on");
Leg.NumColumns = 2;
% ptitle(14,'Daily new cases')
ylim([0 120000])

[Ax(end+1),tidx] = nt(tidx,3);
Sh = plot_mean_var(T.Date,T.ImLossRate,T.ImLossRateStd);
Pl(end+1) = Sh(1);
Pl(end+1) = plot(T.Date,T.ImGainRate,'Color',Color_2);
% Pl(end+1) = plot(T0.Date,T0.ImLossRate,'--','Color',Pl(end-1).Color);
% Pl(end+1) = plot(T0.Date,T0.ImGainRate,'--','Color',Pl(end-1).Color);
Leg(end+1) = legend([Pl(end-1),Sh(4),Pl(end)],{ ...
    'Rate of immunity loss ($\omega$),', ...
    'end its uncertainty ($\pm 95\%$ CI)', ...
    'Immunization rate ($\nu$)'}, ...
    "Interpreter","latex","FontSize",args.FontSize,"Location","northwest");
% ptitle(12,'Estimated waning (proportion of \textbf{R} loosing immunity within 24 hours) and immunization (proportion of \textbf{S} gaining immunity within 24 hours)')
ylim([0,0.02])

T.Rc(end) = T.Rc(end-1);
[Ax(end+1),tidx] = nt(tidx,3);
Pl(end+1) = plot(T.Date,T.Rt,'Color',[0.8500 0.3250 0.0980]);
% Pl(end+1) = plot(T.Date,T.Rc,'k');
Sh = plot_mean_var(T.Date,T.Rc,T.OutputStd(:,J.Rt),[0 0 0],"PMLineStyle",'-',"LineWidth",1);
Pl(end+1) = Sh(1);
% Pl(end+1) = plot(T0.Date,T0.Rc,'k--');
Pl_Unc = area(datetime(1991,03,20) + [0 1],[1 1],'FaceColor',[0 0 0],'FaceAlpha',0.2,'EdgeColor',[0,0,0]);
Leg(end+1) = legend([Pl(end-1:end),Pl_Unc], {...
    '$R_{\rm{c}}$ (OWID)', ...
    '$R_{\rm{c}}$ (computed)', ...
    'and its uncertainty ($\pm 95\%$ CI)'
    ... '$R_{\rm{c}}$ (without)', ...
    }, "Interpreter","latex","FontSize",args.FontSize,"Location","northeast", ...
    "Box","on");
Leg(end).NumColumns = 3;
ylim([0,2.5])

[Ax(end+1),tidx] = nt(tidx,2);
Pl(end+1) = plot(T.Date,1-T.Stringency);
% legend('Stringency idx.', ...
%     "Interpreter","latex","FontSize",args.FontSize,"Location","northwest", ...
%     "Box","off")
text(XLim(1),0.2,'~~Registered stringency index~~',"Interpreter","latex","FontSize",args.FontSize,"HorizontalAlignment","left")
Ax(end).XTickLabel = [];
ylim([0 1])

% Transmission rate
[Ax(end+1),tidx] = nt(tidx,3);
% Sh = plot_interval(T.Date,T.TrRate_bounds(:,1),T.TrRate_bounds(:,2),[0.8500 0.3250 0.0980],"FaceAlpha",0.1);
Pl(end+1) = plot(T.Date,T.TrRate_Ref,'Color',[0.8500 0.3250 0.0980]);
% Pl(end+1) = plot(T.Date,T.TrRate,'k');
Sh = plot_mean_var(T.Date,T.TrRate,T.TrRateStd,[0 0 0],"PMLineStyle",'-',"LineWidth",1);
Pl(end+1) = Sh(1);
% Pl(end+1) = plot(T0.Date,T0.TrRate,'k--');
% Pl(end+1) = plot(T.Date,beta_Is,'Color',[0.9290 0.6940 0.1250]);
Pl_Unc = area(datetime(1991,03,20) + [0 1],[1 1],'FaceColor',[0 0 0],'FaceAlpha',0.2,'EdgeColor',[0,0,0]);
Leg(end+1) = legend([Pl(end) Pl_Unc Pl(end-1)], { ... Sh(3)
    'Transmission rate $\beta$, '
    'its uncertainty ($\pm 95\%$ CI)'
    'and its initial estimate'
    ... 'Admissible range for $\beta$'
    }, ...
    "Interpreter","latex","FontSize",args.FontSize,"Location","northwest", ...
    "NumColumns",3,"Box","on");
% ptitle(14,'Transmission rate of the pathogen')
ylim([0,2.5])

% Hospitalized
[Ax(end+1),tidx] = nt(tidx,3);
Pl(end+1) = plot(T.Date,T.H_off_ma,'Color',[0.8500 0.3250 0.0980]);
% Pl(end+1) = plot(T.Date,T.H,'k');
Sh = plot_mean_var(T.Date,T.H,real(T.StateStd(:,J.H)),[0 0 0],"PMLineStyle",'-',"LineWidth",1);
Pl(end+1) = Sh(1);
% Pl(end+1) = plot(T0.Date,T0.H,'k--');
Pl_Unc = area(datetime(1991,03,20) + [0 1],[1 1],'FaceColor',[0 0 0],'FaceAlpha',0.2,'EdgeColor',[0,0,0]);
Leg(end+1) = legend([Pl(end-1) Pl(end) Pl_Unc], {
    'Hospital load (official data)'
    'Hospital load (computed)'
    'and its uncertainty ($\pm 95\%$ CI)'
    }, ...
    "Interpreter","latex","FontSize",args.FontSize,"Location","northeast", ...
    "NumColumns",3,"Box","on");
% ptitle(args.FontSize,'Deceased people')
Ax(end).YLim(1) = 0;

% Infectious people
[Ax(end+1),tidx] = nt(tidx,4);
% plot(T.Date,T.Szennyviz_Nyers,'--','Color',[0.4660 0.6740 0.1880]);
Pl(end+1) = plot(T.Date,T.Szennyviz .* T.Szennyviz_Mtp,'Color',[0.8500 0.3250 0.0980]);
% Pl(end+1) = plot(T.Date,T.Szennyviz_Nv .* T.Szennyviz_Mtp);
Pl(end+1) = plot(T.Date,T.Infectious,'k');
plot(T0.Date,T0.Infectious,'k--');
Leg(end+1) = plegend(args.FontSize,"northwest",...
    ... "Gene copy conc. (raw);" + ...
    "Gene copy concentration (scaled);" ...
    ... + LegEntries_SzV + ";" ...
    + "Infectious (with wastewater);Infectious (without wastewater)");
Leg(end).NumColumns = 1;
% ptitle(14,'Number infectious people compared with the pathogen concentration of waste water')

[Ax(end+1),tidx] = nt(tidx,3);
Pl(end+1) = plot(T.Date,T.New_Cases);
Pl(end+1) = plot(T.Date,T.New_ImLoss);
Pl(end+1) = plot(T.Date,T.New_ImGain);
Leg(end+1) = legend( ...
    'Daily new cases', ...
    'Daily new immunity loss', ...
    'Daily new immunization', ...
    "Interpreter","latex","FontSize",args.FontSize,"Location","northwest");

[Ax(end+1),tidx] = nt(tidx,4);
Pl(end+1) = plot(T.Date,T.Cum_Cases);
Pl(end+1) = plot(T.Date,T.Cum_ImLoss);
Pl(end+1) = plot(T.Date,T.Cum_ImGain);
Pl(end+1) = plot(T.Date,T.S);
Pl(end+1) = plot(T.Date,T.R);
ylim([0 30000000])
Leg(end+1) = legend( ...
    'Cumulative cases', ...
    'Cumulativ immunity loss', ...
    'Cumulativ immunization', ...
    'Susceptible population', ...
    'Protected population', ...
    "Interpreter","latex","FontSize",args.FontSize,"Location","northwest", ...
    "NumColumns",3);

% -------

Link1 = linkprop(Pl,"LineWidth");
Link2 = linkprop(Ax,"XLim");
Pl(1).LineWidth = 1.2;
Ax(1).XLim = XLim;
drawnow 

return

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
