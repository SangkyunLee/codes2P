function Params = loadVStimParams(Params, pdt)
% function Params = loadVStimParams(Params, pdt)
% INPUT:
% Params collected and pdt: photodiode object
% OUTPUT:
% Store VStimparam structure within Params
% 01-21-2018 Sangkyun Lee

% get the VStim file path and name
files = Params.files;
mainpath = files.mainpath;
if isfield(files,'stim_subpath')
    stim_subpath= files.stim_subpath;
else
    stim_subpath =[];
end
stim_fn1 = files.stim_fn1;
if isfield(files,'stim_mainpath')
    stim_mainpath= files.stim_mainpath;
    stim_fn = fullfile(stim_mainpath,stim_subpath,stim_fn1);
else
    stim_fn = fullfile(mainpath,stim_subpath,stim_fn1);
end




%% loading stimulus parameters

load(stim_fn)
try
    expType = stim.params.constants.expType;
    stimparam = feval(['load_stimparam_' expType],stim.params,Params);
catch 
    error('no stim structure loaded');
end

stim_samp = stimparam.stim_samplesinNI;
if isfield(stimparam, 'blank_samplesinNI')
    blank_samp = stimparam.blank_samplesinNI;
else
    blank_samp =0;
end



%% identify stimulus event from photodiode

t = get(pdt,'t'); % get time of photodiode
y = get(pdt,'y'); % get photodiode  signal

sampf = Params.samplingfreq_NI;
% min interval between two photodiode-event
if blank_samp>0
    minint = min(stim_samp,blank_samp)/sampf; 
else 
    minint = stim_samp/sampf;
end
[tinxe,evtime]=pdt.get_tstmp(minint,[]);

figure; 
t1 = 40;

subplot(2,1,1);
inx1 = t<t1;
inx2 =find(evtime<t1);
plot(t(inx1),y(inx1)); 
hold on; plot(evtime(inx2),max(y)*ones(size(inx2)),'r.');
subplot(2,1,2);
inx1 = t>t(end)-t1;
inx2 =find(evtime>(t(end)-t1));
plot(t(inx1),y(inx1)); 
hold on; plot(evtime(inx2),max(y)*ones(size(inx2)),'r.');



% get stim on and off time

if blank_samp>0 
    lenstim = floor(length(tinxe)/2);
    stimparam.stim_onset_tinxNI = tinxe(1:2:lenstim*2);
    stimparam.blank_onset_tinxNI = tinxe(2:2:lenstim*2);
else
    stimparam.stim_onset_tinxNI = tinxe;
    stimparam.blank_onset_tinxNI = [];    
end
title(sprintf('No. Identified PDT event:%d', length(stimparam.stim_onset_tinxNI)));
check_valid(stimparam);


Params.VStimparam = stimparam;
end

function check_valid(stimparam)
    nrep = stimparam.repetitions;
    stim_onset_tinxNI = stimparam.stim_onset_tinxNI;
    
    switch stimparam.expType
        case{'GratingExperiment2PhotonbySang','SquareMappingExperiment2Photon'}
            if length(stim_onset_tinxNI)<nrep
                warning('Incomplete session: %d(performed), %d(designed)',...
                    length(stim_onset_tinxNI), nrep);
            elseif length(stim_onset_tinxNI)>nrep
                warning('Overcomplete session: %d(performed), %d(designed)',...
                    length(stim_onset_tinxNI), nrep);
            end

        case{'RFmappingBarDrift'}
            nbarloc = length(stimparam.Barcenters);            
            if stimparam.postStimulusTime>0                
                nevt = nbarloc+1;
            else
                nevt = nbarloc;
            end                
            estnrep = length(stim_onset_tinxNI)/nevt;
            if ~(estnrep == nrep)
                warning('estimated no. sweep(%.1f), designed no. sweep(%d)',...
                    estnrep,nrep);                
            end
            
        otherwise
            error('PDT validation code not implemented:%s',stimparam.expType);
    end

end