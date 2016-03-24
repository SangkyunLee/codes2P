function timeinfo = gen_stimtime_GratingExperiment2PhotonbySang(Params)
% function timeinfo = gen_stimtime_GratingExperiment2PhotonbySang(Params)
% Sangkyun Lee 08-30-2013

nstiminblk = 1
nblks0 = Params.stimparam.repetitions;


frame_start = Params.mirror_frame_start;
stimtime = zeros(size(Params.timeNI));
stimtime(Params.stim_onset_tstampinNI)=Params.stimparam.cond.StimEventLUT(1:length(Params.stim_onset_tstampinNI),end);

stimtimelen = length(stimtime);
frametimelen = length(frame_start);


marker = find(diff(stimtime)>0);
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
        cond = stimtime(stpoint);
        stimtime(stpoint:edpoint) = cond;
    end
end

tlen=max(length(frame_start),length(stimtime));
timeinfo.frame_start=zeros(tlen,1);
timeinfo.stimtime=zeros(tlen,1);
timeinfo.frame_start(1:length(frame_start))=frame_start(:);
timeinfo.stimtime(1:length(stimtime))=stimtime(:);
