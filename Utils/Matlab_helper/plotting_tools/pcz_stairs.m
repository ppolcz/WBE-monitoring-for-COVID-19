function [Ph,Pv,Pd] = pcz_stairs(x,y,varargin)
%%
%  File: pcz_stairs.m
%  Directory: utilities/plotting_tools
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2020. May 01. (2019b)
% 

args.Color = pcz_get_plot_colors([],1);

args.HLineStyle = '-';
args.HLineWidth = 2;

args.VLineStyle = ':';
args.VLineWidth = 1;

args.MarkerStyle = '.';
args.MarkerSize = 10;

args = parsepropval(args,varargin{:});

%%

% fig = figure(1);
% delete(get(fig,'Children'));
% ax = axes('Parent',fig);

h = sum(diff(x)/(numel(x)-1));

X = [x(:) ; x(end)+h];
Y = [y(:) ; nan];

XX = [X(:)' ; X(:)' ; X(:)'];
XX = XX(3:end-2)';

YY = [ Y(:)' ; Y(:)' ; nan*Y(:)'];
YY = YY(1:end-4)';

Ph = plot(XX,YY,args.HLineStyle,'Color',args.Color,'LineWidth',args.HLineWidth);
hold on

x = x(:);
XX = [x(2:end)' ; x(2:end)' ; nan*x(2:end)'];
XX = XX(1:end-1)';

YY = [ y(:)' ; nan*y(:)' ; y(:)'];
YY = YY(3:end-2)';

Pv = plot(XX,YY,args.VLineStyle,'Color',args.Color,'LineWidth',args.VLineWidth);

Pd = plot(x,y,args.MarkerStyle,'Color',args.Color,'MarkerSize',args.MarkerSize);

end