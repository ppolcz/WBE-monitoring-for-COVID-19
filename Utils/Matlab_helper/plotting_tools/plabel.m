function [ret] = plabel(AxisNr,FontSize,str,varargin)
%%
%  File: plabel.m
%  Directory: /home/ppolcz/T/_Epid/Utils/Matlab_helper/plotting_tools
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. July 11. (2022a)
%

if ~isempty(varargin) && iscell(varargin{1})
    args = varargin{1};
    varargin = varargin(2:end);
    str = sprintf(str,args{:});
end

args = {str,"FontSize",FontSize,"Interpreter","latex"};
args = [ args varargin ];

switch AxisNr
    case 1
        xlabel(args{:});
    case 2
        ylabel(args{:});
    case 3
        zlabel(args{:});
end

end