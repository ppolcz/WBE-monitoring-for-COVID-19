
addpath(genpath('Utils'))

CasADi_Folder = "Utils/CasADi-Matlab";
if ~exist(CasADi_Folder,"dir")

    GitGub_URL = "https://github.com/casadi/casadi/releases/download/3.6.1/casadi-3.6.1-linux64-matlab2018b.zip";
    ZIP = "casadi.zip";
    
    fprintf('Downloading casadi-3.6.1-linux64-matlab2018b.zip...\n')
    websave(ZIP,GitGub_URL);

    fprintf('Unpacking casadi-3.6.1-linux64-matlab2018b.zip...\n')
    unzip(ZIP,CasADi_Folder);

    fprintf('Deleting casadi-3.6.1-linux64-matlab2018b.zip...\n')
    delete(ZIP);

end
fprintf('Addig CasADi to the Matlab path...\n')
addpath(genpath(CasADi_Folder))
