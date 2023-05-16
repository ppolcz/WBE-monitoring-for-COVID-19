function G_reset(flags)
%%
%  File: G_reset.m
%  Directory: 7_ftools/utilities/global
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2020. April 03. (2019b)
%  Major review on 2020. May 19. (2019b)
%
%  Flags: from least significant digit to most significant:
%  E.g. flag = 01101
%                 ││
%                 │└─ verbosity (0:false, 1:true, 2:do not change)
%                 └── scope depth (0:reset, 1:reset if not set)
%

global SCOPE_DEPTH LATEX_EQNR

if nargin < 1
    flags = 1;
end

verbosity = get_digit(flags,1);
ScopeDepth_resetif = get_digit(flags,2);

G_VERBOSE(verbosity);

if ~ScopeDepth_resetif || isempty(SCOPE_DEPTH)
    SCOPE_DEPTH = 0;
end

if isempty(SCOPE_DEPTH) || ScopeDepth_resetif
    SCOPE_DEPTH = 0;
end

LATEX_EQNR = 0;

end

function digit = get_digit(number,pos)

    digit = ( mod(number,10^pos) - mod(number,10^(pos-1)) ) / 10^(pos-1);
end
