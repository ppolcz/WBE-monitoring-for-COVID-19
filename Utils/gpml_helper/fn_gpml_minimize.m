function ret = fn_gpml_minimize
%%
%  File: fn_gpml_minimize.m
%  Directory: 5_Sztaki20_Main/Utils/gpml_helper
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. February 02. (2020b)
%

persistent gpml_minimize

if isempty(gpml_minimize)
    
    list = which('-all','minimize');
    dir = list(contains(list,'gpml'));
    
    if ~isempty(dir)
        dir = fileparts(dir{1});
        act_dir = pwd;
        cd(dir)
        gpml_minimize = str2func('minimize');
        cd(act_dir)
    else        
        error 'GPML function `minimize` not found'
    end
        
end

ret = gpml_minimize;

