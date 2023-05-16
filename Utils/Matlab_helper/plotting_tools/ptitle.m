function [ret] = ptitle(FontSize,str,varargin)
%%
%  File: ptitle_tex.m
%  Directory: /home/ppolcz/T/_Epid/Utils/Matlab_helper/plotting_tools
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. July 11. (2022a)
%
% Examples:
%  ptitle_simple(14,'Computed infected compartments and new cases')
%  ptitle_simple(14,'Computed infected \n compartments and new cases')
%  ptitle_simple(14,'Computed infected \\n %s compartments and new cases',{'asd'})


if ~isempty(varargin) && iscell(varargin{1})
    args = varargin{1};
    varargin = varargin(2:end);
    str = sprintf(str,args{:});
end

entries = strsplit(str,'\s*\\n\s*','DelimiterType','RegularExpression');

title(entries,"FontSize",FontSize,'Interpreter','latex')

end