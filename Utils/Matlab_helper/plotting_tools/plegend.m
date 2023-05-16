function varargout = plegend(FontSize,Location,str,varargin)
%% 
%  File: plegend_simple.m
%  Directory: [Relative directory in workspace `_` not detected.]
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2022. July 11. (2022a)
% 
% Example:
% 
%  plegend_simple(12,'northwest','$\mathbf{A}$; $\mathbf{B}$ ; $\mathbf{C}$')
% 
%  plegend_simple(12,'west','$%s$; $%s$ ; $%s$',{'\textbf{S}','\textbf{R}','\textbf{H}'})

%% 

if ~isempty(varargin) && iscell(varargin{1})
    args = varargin{1};
    varargin = varargin(2:end);
    str = sprintf(str,args{:});
end

entries = strsplit(str,'\s*;\s*','DelimiterType','RegularExpression');

h = legend(entries{:},varargin{:});
h.Interpreter = 'latex';
h.FontSize = FontSize;
h.Location = Location;

if nargout > 0
    varargout{1} = h;
end

end