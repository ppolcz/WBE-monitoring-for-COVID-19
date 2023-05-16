function varargout = plegend_REGI(varargin)
%% 
%  
%  file:   plegend.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.01.28. Thursday, 19:20:32
%
%% 

h = legend(varargin{:});

h.Interpreter = 'latex';

if nargout > 0
    varargout{1} = h;
end

end