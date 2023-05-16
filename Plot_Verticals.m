function Pl = Plot_Verticals(ax,ticks,varargin)
%%
%  File: plot_verticals.m
%  Directory: /home/ppolcz/T/_Epid/Approach_UIO_then_Opt/Ver_2022_05_31_UIO_then_Opt
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. July 02. (2022a)
%

YLim = ax.YLim;
ticks = ticks(:);

n = numel(ticks);

X = repmat(ticks,[1 3]);
X = X';
X = X(1:end-1);

Y = [ YLim' ; NaN ] .* ones(1,n);
Y = Y(1:end-1);

Pl = plot(ax,X,Y,varargin{:});
ax.YLim = YLim;

end
