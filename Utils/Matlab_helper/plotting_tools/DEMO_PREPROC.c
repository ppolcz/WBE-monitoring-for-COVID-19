%%
%  File: DEMO_PREPROC.m
%  Directory: /home/ppolcz/T/_Epid/Utils/Matlab_helper/plotting_tools
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. July 08. (2022a)
%

%{

#define DEBUG fprintf('Line: %d\n',__LINE__);

DEBUG

%}

G_reset


try c = evalin('caller','persist'); catch; c = []; end
persist = Persist(mfilename('fullpath'), c); clear c; 
persist.backup();
%clear persist

%%
persist.stoplog;
