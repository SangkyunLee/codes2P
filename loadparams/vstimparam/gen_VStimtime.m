function timeinfo= gen_VStimtime(Params, Cond)
% function timeinfo= gen_VStimtime(Params, Cond)
% timeinfo.frame_start:  scanframe time-stamp
% timeinfo.stimtime: stimulus time stamp
% get visual stimulus time 
% 01-21-2018 Sangkyun Lee

stimparam = Params.VStimparam;
stimtime = zeros(size(Params.timeNI));
frame_start = zeros(size(Params.timeNI));




frame_start(Params.tinxf_NI)=10;% oldverion-->Params.mirror_frame_start;


stimstart = stimparam.stim_onset_tinxNI;
% when each trial consists of stim + blank
if isfield(stimparam,'blank_onset_tinxNI') || ~isempty(stimparam.blank_onset_tinxNI)
    stimend = stimparam.blank_onset_tinxNI;
    
else % when each trial consists only of stim
    stimend = stimstart(2:end);
    est_stimend = stimend(end) + round(mean(diff(stimstart)));
    if est_stimend<=length(stimtime)
        stimend(end+1) = est_stimend;
    else
        stimend(end+1) = length(stimtime);
        warning('The last trial(%d) is cut out',length(stimstart));
    end
end

nrep = stimparam.repetitions; % designed repetition
nrep1 = length(stimstart); %executed repettion
stimtime(stimparam.stim_onset_tinxNI)=Cond.condseq(1:nrep1);
if nrep1> nrep    
    fprintf('%d(designed), %d(executed): %d used for data collection',...
        nrep, nrep1, nrep);
    nrep1 = nrep;
end
    
    
for iblk = 1: nrep1
    stpoint = stimstart(iblk);          
    edpoint = stimend(iblk);     
    cond = stimtime(stpoint);
    inxes = stpoint : edpoint;
    stimtime(inxes) = cond;    
end    


tlen=max(length(frame_start),length(stimtime));
timeinfo.frame_start=zeros(tlen,1);
timeinfo.stimtime=zeros(tlen,1);
timeinfo.frame_start(1:length(frame_start))=frame_start(:);
timeinfo.stimtime(1:length(stimtime))=stimtime(:);