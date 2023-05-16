function subl(varargin)
%% 
%  
%  file:   subl.m
%  author: Polcz Péter <ppolcz@gmail.com> 
% 
%  Created on Fri Sep 19 18:04:08 CEST 2014
%
%%

if nargin > 0 && ~isempty(varargin{1})
	fn = varargin{1};
else
    active = matlab.desktop.editor.getActive();
    line = active.Selection(1:2);
    fn = [active.Filename ':' num2str(line(1)) ':' num2str(line(2)) ];
end

args = '';
if (nargin > 1)
    args = strjoin(varargin(2:end));
end

system(['subl ', args, ' ', fn]);

end
