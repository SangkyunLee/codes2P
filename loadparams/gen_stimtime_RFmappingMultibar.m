function timeinfo = gen_stimtime_RFmappingMultibar(Params)
% function timeinfo = gen_stimtime_RFmappingMultibar(Params)
% output:
%    timeinfo.frame_start
%    timeinfo.stimtime1: image order in the time line of the NI card
%    timeinfo.stimtime2: Orientation order in the time line of the NI card
% Sangkyun Lee 10-15-2013


nblks0 = Params.stimparam.nblks;


frame_start = Params.mirror_frame_start;
stimtime1 = zeros(size(Params.timeNI));
stimtime1(Params.stim_onset_tstampinNI) = Params.stimparam.imgOrder(1:length(Params.stim_onset_tstampinNI));
stimtime2 = zeros(size(Params.timeNI));
stimtime2(Params.stim_onset_tstampinNI) = Params.stimparam.ORIOrder(1:length(Params.stim_onset_tstampinNI));

stimtimelen = length(stimtime1);
frametimelen = length(frame_start);

nstiminblk = Params.stimparam.nstiminblk; %designed No. stim in a block
stimleninblk = nstiminblk; %No. stim run in a block


marker = find(diff(stimtime1)>0);
nblks = length(marker)/nstiminblk;
if nblks < nblks0
    bincomplete =1;
    nstimlastblk = mod(nblks*nstiminblk,nstiminblk);
    nblks =ceil(nblks);
else
    bincomplete =0;
end
for iblk=1:nblks
    if iblk == nblks & bincomplete,
        stimleninblk = nstimlastblk;
    end            
    for istim=1:nstiminblk
        stpoint = marker((iblk-1)*nstiminblk+istim)+1;
        if istim==nstiminblk,
            if bincomplete,
                edpoint = stpoint + 1;
            else
                edpoint = Params.blank_onset_tstampinNI(iblk);
            end
        else                
            edpoint = marker((iblk-1)*nstiminblk+istim+1);                
        end
        cond = stimtime1(stpoint);
        stimtime1(stpoint:edpoint) = cond;
        cond = stimtime2(stpoint);
        stimtime2(stpoint:edpoint) = cond;
    end
end

tlen=max(length(frame_start),length(stimtime1));
timeinfo.frame_start=zeros(tlen,1);
timeinfo.stimtime1=zeros(tlen,1);
timeinfo.stimtime2=zeros(tlen,1);
timeinfo.frame_start(1:length(frame_start))=frame_start(:);
timeinfo.stimtime1(1:length(stimtime1))=stimtime1(:);
timeinfo.stimtime2(1:length(stimtime2))=stimtime2(:);
