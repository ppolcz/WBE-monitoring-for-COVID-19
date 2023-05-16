function [getvalue] = G_VERBOSE(setvalue)
%% G_VERBOSE
%  
%  File: G_VERBOSE.m
%  Directory: 7_ftools/utilities/global
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. April 03. (2019b)
%  Major review on 2020. May 19. (2019b)
%

global VERBOSE

if isempty(VERBOSE)
    VERBOSE = true;
end

if nargin == 1
    getvalue = VERBOSE;
    if setvalue < 2
        VERBOSE = logical(setvalue);
    end
    return
elseif ~isscalar(VERBOSE)
    VERBOSE = false;
end

getvalue = VERBOSE;

end