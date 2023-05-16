% Plot colors (X11 color names)

Color_1 = [0 0.4470 0.7410];
Color_2 = [0.8500 0.3250 0.0980];
Color_3 = [0.9290 0.6940 0.1250];
Color_4 = [0.4940 0.1840 0.5560];
Color_5 = [0.4660 0.6740 0.1880];
Color_6 = [0.3010 0.7450 0.9330];
Color_7 = [0.6350 0.0780 0.1840];

Color_Black           = pcz_hex2rgb("#000000");
Color_Light_Black     = pcz_hex2rgb("#808080");
Color_Blue            = pcz_hex2rgb("#0000FF");
Color_Dark_Blue       = pcz_hex2rgb("#00008B");
Color_Light_Blue      = pcz_hex2rgb("#ADD8E6");
Color_Brown           = pcz_hex2rgb("#A52A2A");
Color_Dark_Brown      = pcz_hex2rgb("#5C4033");
Color_Light_Brown     = pcz_hex2rgb("#996600");
Color_Buff            = pcz_hex2rgb("#F0DC82");
Color_Dark_Buff       = pcz_hex2rgb("#976638");
Color_Light_Buff      = pcz_hex2rgb("#ECD9B0");
Color_Cyan            = pcz_hex2rgb("#00FFFF");
Color_Dark_Cyan       = pcz_hex2rgb("#008B8B");
Color_Light_Cyan      = pcz_hex2rgb("#E0FFFF");
Color_Gold            = pcz_hex2rgb("#FFD700");
Color_Dark_Gold       = pcz_hex2rgb("#EEBC1D");
Color_Light_Gold      = pcz_hex2rgb("#F1E5AC");
Color_Goldenrod       = pcz_hex2rgb("#DAA520");
Color_Dark_Goldenrod  = pcz_hex2rgb("#B8860B");
Color_Light_Goldenrod = pcz_hex2rgb("#FFEC8B");
Color_Gray            = pcz_hex2rgb("#808080");
Color_Dark_Gray       = pcz_hex2rgb("#404040");
Color_Light_Gray      = pcz_hex2rgb("#D3D3D3");
Color_Green           = pcz_hex2rgb("#008000");
Color_Dark_Green      = pcz_hex2rgb("#006400");
Color_Light_Green     = pcz_hex2rgb("#90EE90");
Color_Ivory           = pcz_hex2rgb("#FFFFF0");
Color_Dark_Ivory      = pcz_hex2rgb("#F2E58F");
Color_Light_Ivory     = pcz_hex2rgb("#FFF8C9");
Color_Magenta         = pcz_hex2rgb("#FF00FF");
Color_Dark_Magenta    = pcz_hex2rgb("#8B008B");
Color_Light_Magenta   = pcz_hex2rgb("#FF77FF");
Color_Mustard         = pcz_hex2rgb("#FFDB58");
Color_Dark_Mustard    = pcz_hex2rgb("#7C7C40");
Color_Light_Mustard   = pcz_hex2rgb("#EEDD62");
Color_Orange          = pcz_hex2rgb("#FFA500");
Color_Dark_Orange     = pcz_hex2rgb("#FF8C00");
Color_Light_Orange    = pcz_hex2rgb("#D9A465");
Color_Pink            = pcz_hex2rgb("#FFC0CB");
Color_Dark_Pink       = pcz_hex2rgb("#E75480");
Color_Light_Pink      = pcz_hex2rgb("#FFB6C1");
Color_Red             = pcz_hex2rgb("#FF0000");
Color_Dark_Red        = pcz_hex2rgb("#8B0000");
Color_Light_Red       = pcz_hex2rgb("#FF3333");
Color_Silver          = pcz_hex2rgb("#C0C0C0");
Color_Dark_Silver     = pcz_hex2rgb("#AFAFAF");
Color_Light_Silver    = pcz_hex2rgb("#E1E1E1");
Color_Turquoise       = pcz_hex2rgb("#30D5C8");
Color_Dark_Turquoise  = pcz_hex2rgb("#00CED1");
Color_Light_Turquoise = pcz_hex2rgb("#AFE4DE");
Color_Violet          = pcz_hex2rgb("#EE82EE");
Color_Dark_Violet     = pcz_hex2rgb("#9400D3");
Color_Light_Violet    = pcz_hex2rgb("#7A5299");
Color_White           = pcz_hex2rgb("#FFFFFF");
Color_Yellow          = pcz_hex2rgb("#FFFF00");
Color_Dark_Yellow     = pcz_hex2rgb("#FFCC00");
Color_Light_Yellow    = pcz_hex2rgb("#FFFFE0");

New_Cases_Colors = {
    Color_Light_Blue                              "Wild"
    Color_Goldenrod                               "Alpha"
    Color_5*1.2                                   "Delta"
    pcz_colormix(Color_2,0.7)                     "BA.1"
    Color_Silver                                  "BA.2"
    pcz_colormix(Color_Dark_Ivory,0.16,Color_Gold) "BA.5"
    Color_Light_Pink                              "BA.5--Tp1"
    Color_Turquoise                               "BQ.1"
    };

%%

colororder = ([
    Color_1
    Color_2
    Color_3
    Color_4
    Color_5
    Color_6
    Color_7
    Color_Black
    Color_Light_Black
    Color_Blue
    Color_Dark_Blue
    Color_Light_Blue
    Color_Brown
    Color_Dark_Brown
    Color_Light_Brown
    Color_Buff
    Color_Dark_Buff
    Color_Light_Buff
    Color_Cyan
    Color_Dark_Cyan
    Color_Light_Cyan
    Color_Gold
    Color_Dark_Gold
    Color_Light_Gold
    Color_Goldenrod
    Color_Dark_Goldenrod
    Color_Light_Goldenrod
    Color_Gray
    Color_Dark_Gray
    Color_Light_Gray
    Color_Green
    Color_Dark_Green
    Color_Light_Green
    Color_Ivory
    Color_Dark_Ivory
    Color_Light_Ivory
    Color_Magenta
    Color_Dark_Magenta
    Color_Light_Magenta
    Color_Mustard
    Color_Dark_Mustard
    Color_Light_Mustard
    Color_Orange
    Color_Dark_Orange
    Color_Light_Orange
    Color_Pink
    Color_Dark_Pink
    Color_Light_Pink
    Color_Red
    Color_Dark_Red
    Color_Light_Red
    Color_Silver
    Color_Dark_Silver
    Color_Light_Silver
    Color_Turquoise
    Color_Dark_Turquoise
    Color_Light_Turquoise
    Color_Violet
    Color_Dark_Violet
    Color_Light_Violet
    Color_White
    Color_Yellow
    Color_Dark_Yellow
    Color_Light_Yellow
    ]);

    Colors = colororder;

%%

function plotColors
%%
    Plot_Colors

    fig = figure(41231);
    fig.Position(3) = 1960;
    Tl = tiledlayout(1,1,'Padding','compact');
    nexttile, hold on;

    Colors = colororder;

    for r=1:height(Colors)
        x = linspace(0,r,500);
        y = sqrt(r.^2-x.^2);
        plot(x,y,'LineWidth',15,'Color',Colors(r,:))
    end

    grid on
    ylim([0,1])
    xlim([0,height(Colors)])
    xticks(1:height(Colors))

end
