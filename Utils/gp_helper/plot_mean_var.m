function Sh = plot_mean_var(tt,xx,xx_std,FaceColor,args)
arguments
    tt,xx,xx_std,
    FaceColor = [0 0.4470 0.7410]
    args.LineColor = FaceColor;
    args.LineWidth = 1.2;
    args.LineStyle = '-';
    args.PMLineStyle = 'none';
    args.PlotMean = true;
    args.Alpha = 2;
    args.FaceAlpha = 0.2;
end

    if args.PlotMean
        Pl_Mu = plot(tt,xx,'Color',args.LineColor,'LineStyle',args.LineStyle,"LineWidth",args.LineWidth);
    end

    Sh = shade(tt,xx'+args.Alpha*xx_std',tt,xx-args.Alpha*xx_std, ...
        'FillType',[1 2;2 1], ...
        'FillColor',FaceColor);

    try
        Sh(3).FaceAlpha = args.FaceAlpha;
        Sh(1).LineStyle = args.PMLineStyle;
        Sh(1).Color = args.LineColor;
        Sh(2).LineStyle = args.PMLineStyle;
        Sh(2).Color = args.LineColor;
    end

    if args.PlotMean
        Sh = [Pl_Mu ; Sh];
    end

end