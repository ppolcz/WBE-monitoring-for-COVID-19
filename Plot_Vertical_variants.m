function [ret] = Plot_Vertical_variants(ax,Q,args)
arguments
    ax,Q
    args.Except = []
end
%%
%  File: plot_vertical_variants.m
%  Directory: /home/ppolcz/T/_Epid/RecPred_UIO_then_Opt/Ver_2022_07_12_UIO_then_Opt
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. July 15. (2022a)
%

% 2022.09.22. (szeptember 22, csütörtök), 13:07
Variants = Q.Properties.RowNames;
Variants = strrep(Variants,"Original","Wild");
Variants = strrep(Variants,"OmicronBA","BA.");
Variants = strrep(Variants,"OmicronBQ","BQ.");
Variants = strrep(Variants,"BA.5_Tp1","BA.5\,$+$\,School");

% 2022.09.22. (szeptember 22, csütörtök), 13:21
Idx = find(~ismember(Variants,["Transient","Future" , args.Except]))';
for i = Idx
    if ax.XLim(1) <= Q.Date(i) && Q.Date(i) <= ax.XLim(2)

        % 2022.09.22. (szeptember 22, csütörtök), 13:07
        % Text = sprintf('%s (%s)',Variants(i),datestr(Q.Date(i),'yyyy-mmm-dd'));
        Text = sprintf('%s',Variants(i));

        Plot_Vertical_milestone(ax,Q.Date(i),Text,'k');

        if ~isempty(ax.Legend)
            Idx = strcmp("data1",ax.Legend.String);
            ax.Legend.String(Idx) = [];
        end
    end
    
end

end