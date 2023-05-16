function [Colors,Axh,Leg] = pcz_prelegend(Leg_Entries,indices,varargin)
%%
%  File: pcz_prelegend.m
%  Directory: 7_ftools/ftools/v12/utilities/plotting_tools
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. June 18. (2020a)
%

N = numel(Leg_Entries);
Colors = num2cell(pcz_get_plot_colors([],1:N),2);

hold off
Axh(1) = plot(0,0,'Color',Colors{indices(1)},varargin{:});
hold on;
for i = 2:numel(indices)
    Axh(i) = plot(0,0,'Color',Colors{indices(i)},varargin{:});
end

Leg = plegend(Leg_Entries{:});

end