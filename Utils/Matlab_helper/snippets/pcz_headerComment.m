function [varargout] = pcz_headerComment(varargin)
%%
%  file:   pcz_headerComment.m
%  author: Polcz Peter (ppolcz@gmail.com)
%
%  Created on 2016.01.31. Sunday, 12:32:05
%  Reviewed on 2017. October 18. [to mlx]
%  Minor review on 2020. May 19. (2019b)
% 

editor = 0;
fn = '';
dir = '';
if nargin > 0 && ischar(varargin{1})
    fn = varargin{1};
    [reldir,~,~] = fileparts(fn);
    dir = strrep(pwd,proot,'');
    if ~isempty(reldir)
        dir = [ dir '/' reldir ];
    end
elseif nargin > 0 && isstruct(varargin{1})
    f = varargin{1};
    fn = f.fn;
    dir = f.reldir;
elseif nargin == 0
    active = matlab.desktop.editor.getActive;
    f = pcz_mfilename(active.Filename);
    fn = f.fname;
    dir = f.reldir;
    editor = 1;
end


ret = sprintf([...
    '%%%%\n'...
    ... '%%  File: ' fn '\n'...
    ... '%%  Directory: ' dir '\n'...
    '%%  Author: Peter Polcz (ppolcz@gmail.com) \n'...
    ... '%%  \n'...
    '%%  Created on ' pcz_fancyDate('informative') ' (' version('-release') ')\n'...
    '%%\n\n']);

if editor
    ret_text = sprintf([...
        ... ' File: ' fn '\n'...
        ... ' Directory: ' dir '\n'...
        ' Author: Peter Polcz (ppolcz@gmail.com) \n\n'...
        ' Created on ' pcz_fancyDate('informative') ' (' version('-release') ')\n']);
else
    ret_text = sprintf([...
        ... 'File: ' fn '\n'...
        ... 'Directory: ' dir '\n'...
        'Author: Peter Polcz (ppolcz@gmail.com) \n\n'...
        'Created on ' pcz_fancyDate('informative') ' (' version('-release') ')\n']);
end
    
if nargout > 0
    varargout{1} = ret;
else
    clipboard('copy', ret_text)
end

end