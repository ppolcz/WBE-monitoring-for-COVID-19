%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2023. April 26. (2022b)
function Plot_WData_City(R0,R1,T_by_city,city,args)
arguments
    R0,R1
    T_by_city
    city = 'szeged'
    args.XLim = [ datetime(2020,08,31) , datetime(2023,01,01) ]
    args.FontSize = 12
end

Plot_Colors

XLim = args.XLim;
FontSize = args.FontSize;
Date = XLim(1):XLim(2);

City = [ upper(city(1)) city(2:end) ];

fig = figure(734);
fig.Position(3:4) = [619 447];
Tl = tiledlayout(3,1,"TileSpacing","tight","Padding","compact");

ldx0 = R0.Varos == city & XLim(1) <= R0.Date & R0.Date < XLim(2);
ldx1 = R1.Varos == city & XLim(1) <= R1.Date & R1.Date < XLim(2);

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

plot(XLim,[1 1]*median(R0.LoD),'Color',Color_5,'LineWidth',2,...
    'DisplayName','Limit of detection')

box on
xlim(XLim)
ylim([10e2 10e6])
ylabel('Genome copy concentration [GC/L]',"Interpreter","latex","FontSize",args.FontSize)
Leg = legend('Location','northwest',"Interpreter","latex","FontSize",args.FontSize);

Ax.XTickLabel = [];

drawnow

Ax(2) = nexttile; hold on, grid on
% Xl = xline(Date(day(Date) == 1),'Color',[1 1 1]*0.75);
% Ax(2).YGrid = "on";
% for xl = Xl'
%     xl.HandleVisibility = "off";
% end

Cp = T_by_city.(city);

ldx = ~isnan(Cp);

plot(R0.Date(ldx0),R0.Cp(ldx0),'.-','MarkerSize',20,'Color',[1 1 1]*0.7,'HandleVisibility','off')
plot(R1.Date(ldx1),R1.Cp(ldx1),'.-','MarkerSize',20,'Color',Color_1,'HandleVisibility','off')
plot(T_by_city.Date,T_by_city.(city),'LineWidth',3,'Color',Color_2,'DisplayName','Moving median')
Leg = legend('Location','east',"Interpreter","latex","FontSize",args.FontSize);

box on
xlim(XLim)
ylim([0 40])
ylabel('Cycles ($C_\mathrm{q}$)',"Interpreter","latex","FontSize",args.FontSize)

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
fname = DIR + "/PCR-Results_" + City + "-" + tag + string(datetime(date,'Format','uuuu-MM-dd')) + ".pdf";
% exportgraphics(fig,fname,"ContentType","vector")


end