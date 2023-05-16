function [Pl,Tx] = Plot_Vertical_milestone(ax,xval,Text,varargin)
%%
%  File: plot_verticals.m
%  Directory: /home/ppolcz/T/_Epid/Approach_UIO_then_Opt/Ver_2022_05_31_UIO_then_Opt
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. July 02. (2022a)
%

YLim = ax.YLim;
Pl = plot(ax,[xval xval],YLim,varargin{:});
ax.YLim = YLim;

Vertical = true;
if startsWith(Text,"[HOR]")
    Vertical = false;
    Text = strrep(Text,"[HOR]","");
end

VERTICAL = false;
if startsWith(Text,"[VER]")
    VERTICAL = true;
    Text = strrep(Text,"[VER]","");
end

% 2022.09.22. (szeptember 22, csütörtök), 13:33
if (strlength(Text) < 20 || ~Vertical) && ~VERTICAL
    % Horizontal
    Rotation = 0;
    Text = newline + "~" + Text;
    VerticalAlignment = "top";
    HorizontalAlignment = "left";
else
    % Vertical
    Rotation = -90;
    Text = "~~" + Text;
    VerticalAlignment = "bottom";
    HorizontalAlignment = "left";
end

Tx = text(ax,xval,YLim(2),Text, ...
    "Rotation",Rotation, ...
    "HorizontalAlignment",HorizontalAlignment, ...
    "VerticalAlignment",VerticalAlignment, ...
    "FontSize",10, ...
    "Interpreter","latex");

end
