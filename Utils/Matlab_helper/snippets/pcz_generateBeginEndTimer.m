function varargout = pcz_generateBeginEndTimer
%% 
%  
%  File: pcz_generateBeginEndTimer.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2016.01.17. Sunday, 13:30:51
%  Modified on 2018. April 16.
%

%
%%

% 2021.12.25. (december 25, szombat), 03:18
rng shuffle

var = ['Timer_' pcz_generateString(4, 1) ];
beginning = sprintf('%s = pcz_dispFunctionName;\n', var);
ending = sprintf('\npcz_dispFunctionEnd(%s);\nclear TMP_*\n', var);

%%
if nargout == 1
    varargout{1} = [beginning ending];
elseif nargout == 2
    varargout{1} = beginning;
    varargout{2} = ending;
else
    beginning = sprintf('pcz_dispFunction('''')\n%s = pcz_dispFunctionName;\n', var);
    ending = sprintf('pcz_dispFunctionEnd(%s);', var);
    clipboard('copy', [beginning ending])
end
    
end
