function T = tb_expand_timerange(T,trange,args)
arguments
    T,trange
    args.TimeVariables = ["Date","Day","Day_MPC"]
end
%%
%  File: pf_adjust_tables.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study19_timedep_newV
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2022. May 05. (2022a)
%

%%

% Generate rows before
if trange(1) < T.Date(1)
    Npp = days(T.Date(1) - trange(1));

    T(1+Npp:end+Npp,:) = T;
    T(1:Npp,:) = repmat(T(1,:),[Npp 1]);

    for fn = args.TimeVariables
        if ismember(fn,T.Properties.VariableNames)
            T.(fn)(1:Npp) = T.(fn)(Npp+1) + (-Npp:-1)';
        end
    end

    if isa(T,"timetable")
        T.Properties.RowTimes(1:Npp) = T.Properties.RowTimes(Npp+1) + (-Npp:-1)';
    end
end
    
% Generate rows after
if T.Date(end) < trange(2)
    N = size(T,1);
    Npp = days(trange(2) - T.Date(end));

    T(end+1:end+Npp,:) = repmat(T(end,:),[Npp 1]);
    
    for fn = args.TimeVariables
        if ismember(fn,T.Properties.VariableNames)
            T.(fn)(N+1:end) = T.(fn)(N) + (1:Npp)';
        end
    end

    if isa(T,"timetable")
        T.Properties.RowTimes(N+1:end) = T.Properties.RowTimes(N) + (1:Npp)';
    end
end

end