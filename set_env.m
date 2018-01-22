function parpoolid = set_env(bpar,npool)
if ispc
    addpath(genpath('D:\packages\codes2P'));
    addpath(genpath('Z:\codes2P'));
    addpath(genpath('./helpfun'));

  
elseif isunix
    [~, hn] = system('hostname');
    if strfind(hn,'SLCOMPSS')
        addpath(genpath('/home/slee/data/codes2P'));
        addpath(genpath('./helpfun'));
    else
        addpath(genpath('/home/sl447/codes2P'));
        addpath(genpath('./helpfun'));
    end
end
if bpar && str2double(strtok(version,'.'))>7
    if isunix
        if strfind(hn,'SLCOMPSS')
            if ~exist('npool','var')
                npool=39;       
            end
        else
            if ~exist('npool','var')
                npool=100;       
            end
        end
    else
        npool=6;
    end
else
    
end    
if exist('npool','var')==1
    delete(gcp('nocreate'))
    parpoolid = parpool(npool);
    pause(1);
else
    parpoolid =[];
end
