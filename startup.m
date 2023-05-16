
addpath(genpath('Utils'))

CasADi_Folder = "Utils/CasADi-Matlab";
if ~exist(CasADi_Folder,"dir")

    GitGub_URL = "https://github.com/casadi/casadi/releases/download/3.6.1/casadi-3.6.1-linux64-matlab2018b.zip";
    ZIP = "casadi.zip";
    
    websave(ZIP,GitGub_URL)
    unzip(ZIP,CasADi_Folder)    
    delete(ZIP)

end
addpath(genpath(CasADi_Folder))


GPML_Folder = "Utils/gpml-matlab-master";
if ~exist(GPML_Folder,"dir")

    GitLab_URL = "https://gitlab.com/hnickisch/gpml-matlab/-/archive/master/gpml-matlab-master.zip";
    ZIP = "gpml.zip";
    
    websave(ZIP,GitLab_URL)
    unzip(ZIP,GPML_Folder)    
    delete(ZIP)

end
addpath(genpath(GPML_Folder))
