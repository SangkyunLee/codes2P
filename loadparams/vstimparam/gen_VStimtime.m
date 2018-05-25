function timeinfo= gen_VStimtime(Params, Cond)
% function timeinfo= gen_VStimtime(Params, Cond)
% Cond is struct with 'condseq'
% for condseq, '0' is the blank condition
% timeinfo.frame_start:  scanframe time-stamp
% timeinfo.stimtime: stimulus time stamp
% get visual stimulus time 
% 01-21-2018 Sangkyun Lee

stimparam = Params.VStimparam;
stimtime0 = zeros(size(Params.timeNI));
stimtime = stimtime0;

frame_start = zeros(size(Params.timeNI));
frame_start(Params.tinxf_NI)=10;% oldverion-->Params.mirror_frame_start;


stimstart = stimparam.stim_onset_tinxNI;
nevt1 = length(stimstart); %executed repettion detected by photodiode
stimtime0(stimparam.stim_onset_tinxNI)=Cond.condseq(1:nevt1);

% when each trial consists of stim + blank
if isfield(stimparam,'blank_onset_tinxNI') && ~isempty(stimparam.blank_onset_tinxNI)
    stimend = stimparam.blank_onset_tinxNI;
    
else % when each trial consists only of stim
    
    %Cond.condseq==0 ==> blank condition
    if Cond.condseq(end)>0
        % if the last event is not blank stimulus,
        % approximately estimate the stimevent end
        stimend = stimstart(2:end)-1;
        stimdur = diff(stimstart);
        stimdur = round(mean(stimdur(Cond.condseq(1:end-1)>0)));
        est_stimend = stimend(end) + stimdur;
        if est_stimend<=length(stimtime)
            stimend(end+1) = est_stimend;
        else
            stimend(end+1) = length(stimtime);
            warning('The last trial(%d) is cut out',length(stimstart));
        end
    else
        stimend =[];
    end
end



% designed repetition

nevt = length(Cond.condseq); 

if nevt1> nevt
    fprintf('NEVT:%d(designed), %d(executed): %d used for data collection',...
        nevt, nevt1, nevt);
    nevt1 = nevt;
end

    
    




for iblk = 1: nevt1
    stpoint = stimstart(iblk);     
    if isempty(stimend)
        edpoint = length(stimtime);
    else
        edpoint = stimend(iblk);
    end
    cond = stimtime0(stpoint);
    inxes = stpoint : edpoint;
    stimtime(inxes) = cond;    
end    


tlen=max(length(frame_start),length(stimtime));
timeinfo.frame_start=zeros(tlen,1);
timeinfo.stimtime=zeros(tlen,1);
timeinfo.frame_start(1:length(frame_start))=frame_start(:);
timeinfo.stimtime(1:length(stimtime))=stimtime(:);