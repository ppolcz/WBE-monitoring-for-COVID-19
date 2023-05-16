function [fig,Tl] = Plot_Prediction(T0,T,Q,Qty,Pred,args)
arguments
    T0,T,Q,Qty,Pred
    args.Date_Start_Pred
    args.FontSize = 11;
    args.XLim = [];
    args.TileCols = 1;
    args.ColSpan = 1;
end
%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2022. October 21. (2022b)
%

%% Visualization

Plot_Colors

if isempty(args.XLim)
    args.XLim = [ datetime(2021,11,01) T.Properties.UserData.Date_Last_Szennyviz ];
    % args.XLim = [ datetime(2022,01,01) T.Date(end) ];
end

fig = figure(1312);
fig.Position(3:4) = [574 527];
Tl = tiledlayout(3,args.TileCols,"TileSpacing","tight","Padding","compact","TileIndexing","columnmajor");

Idx0 = T0.Date <= args.Date_Start_Pred;
IdxF = T0.Date >= args.Date_Start_Pred;

LegEntries = {
    'Rate of immunity loss ($\omega_k$)'         'Reconstruction' 'Expectation' 'Possible future scenarios'
    'Transmission rate ($\beta_k$)'             'Reconstruction' 'Expectation' 'Possible future scenarios'
    'Hospitalization probability'   'Reconstruction' 'Expectation' 'Possible future scenarios'
    'New cases ($u_k$)'                     'Reconstruction' 'Prediction' 'Uncertainty'
    'All infected'                  'Reconstruction' 'Prediction' 'Uncertainty'
    'Hospitalized patients'         'Reconstruction' 'Prediction' 'Uncertainty'
    };


tidx = 1;
k = 0;
for q = [1 2 4]
    k = k + 1;
    if q == 1
        Ax = nexttile;
    else
        Ax(end+1) = nexttile;
    end
    hold on;
    grid on;
    box on;

    qty = [Pred{q,:}];
    Mean = mean(qty,2);
    Std = std(qty,0,2);

    fprintf('%s: Mean = %g, Std = %g, Perc = %g, Perc 2sigma = %g\n', ...
        LegEntries{q,1},Mean(end),Std(end),(qty(end,end)-Mean(end))/Mean(end)*100, ...
        round(2*Std(end)/Mean(end)*100))


    plot(datetime(1991,03,20),0,'w','DisplayName','')
    Pl0 = plot(T0.Date,Qty{q,1}(T0),'LineWidth',1.5,'Color',[0 0.4470 0.7410]);
    Pl0.DisplayName = "Reconstruction";

    Xl = xline(args.Date_Start_Pred);
    Xl.HandleVisibility = "on";
    Xl.Color = Color_Light_Red;
    Xl.LineWidth = 2;
    Xl.DisplayName = "Prediction started (vertical line)~~~~~~";
    
    plot(datetime(1991,03,20),0,'w','DisplayName','')    
    Sh = plot_mean_var(T.Date(1:end-1),Mean,Std,[0.9290 0.6940 0.1250],"FaceAlpha",0.2,"LineWidth",1.5);
    Sh(1).DisplayName = "Expectation";
    Sh(2).HandleVisibility = "off";
    Sh(3).HandleVisibility = "off";
    Sh(4).DisplayName = "$95\%$ CI";
    % plot(T.Date(1:end-1),[Dat1{q,:}],'Color',[0.9290 0.6940 0.1250])

    plot(T0.Date,Qty{q,1}(T0),'Color',[0 0.4470 0.7410],"HandleVisibility","off");
    
    title("\textbf{P" + num2str(k) + ".} " + LegEntries{q,1},"FontSize",11,"Interpreter","latex")

    % PlF = plot(T0.Date(IdxF),Qty{q,1}(T0(IdxF,:)),'--','LineWidth',1.5,'Color',[0 0.4470 0.7410]);
    % PlF.DisplayName
    % plot(T0.Date(IdxF),Qty{q,1}(T0(IdxF,:)),'-','Color',[0 0.4470 0.7410]);
%     legend([Pl0 PlF Sh(1) Sh(4)], LegEntries(q,:), ...
%         "Interpreter","latex","FontSize",args.FontSize,"Location","northwest")   

    if q == 1
        % for k = 1:4, plot(datetime(1991,03,20),0,'w','DisplayName','~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'); end
        Leg = legend('Interpreter','latex','Location','northoutside','Box','off','FontSize',12,...
            "NumColumns",2);
    end

    drawnow
end

Link6 = linkprop(Ax,"XLim");
Ax(1).XLim = args.XLim;

% for i = 1:numel(Ax)
%     Pl = plot_vertical_milestone(Ax(i),args.Date_Start_Pred,"Prediction");
%     Pl.HandleVisibility = "off";
% end

% -------

Link2 = linkprop(Ax,"XLim");
Ax(1).XLim = args.XLim;
drawnow 

% return

Date = min(T0.Date(1),T.Date(1)):max(T.Date(end),T0.Date(end));

% -------

for idx = 1:numel(Ax)
    ax = Ax(idx);
    Logger.latexify_axis(ax,args.FontSize);
    % ax.XTick = Date(day(Date) == 1 & mod(month(Date),3) == 1);
    ax.XTick = Date(day(Date) == 1 & mod(month(Date),3) == 0);
    ax.XMinorGrid = 'on';
    ax.XAxis.MinorTick = 'off';
    ax.XAxis.MinorTickValues = Date(day(Date) == 1);    
    drawnow 
end

Today = datetime('today','Format','uuuu-MM-dd');
exportgraphics(fig,"/home/ppolcz/Dropbox/Peti/Munka/01_PPKE_2020/Dokumentaciok/Docs_CsutakB_PhD/07_COVID-Szennyviz/actual/fig/Prediction_v2_" + string(Today) + ".pdf",'ContentType','vector')

end
