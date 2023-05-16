function T = pcz_read_epid_table(UserData,varargin)
    varnames = ["Param","ParamStd","Phases","Waning_Bases"];

    T = readtimetable(varargin{:});

    for vn = varnames
        merge_vns = T.Properties.VariableNames(~cellfun(@isempty,regexp(T.Properties.VariableNames,"^"+vn+"_\d*$")));
        T = mergevars(T,merge_vns,"NewVariableName",vn);
    end
    T.Properties.UserData = UserData;
end
