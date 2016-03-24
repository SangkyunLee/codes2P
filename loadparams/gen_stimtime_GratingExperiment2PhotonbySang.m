function timeinfo = gen_stimtime_GratingExperiment2PhotonbySang(Params)
% % function timeinfo = gen_stimtime_GratingExperiment2PhotonbySang(Params)


stimtime = zeros(size(Params.timeNI));
stimtime(Params.stim_onset_tstampinNI)=Params.stimparam.cond.StimEventLUT(1:length(Params.stim_onset_tstampinNI),end);

stimtimelen=length(stimtime);
frame_start = Params.mirror_frame_start;
frametimelen=length(frame_start);
if ~isfield(Params.stimparam,'nstiminblk')
    nstiminblk =1;
else
    nstiminblk = Params.stimparam.nstiminblk; %designed No. stim in a block
end
stimleninblk = nstiminblk; %No. stim run in a block

marker = find(diff(stimtime)>0);
nblks = length(marker)/nstiminblk;

if isfield(Params.stimparam,'nblks') 
    nblks2 = Params.stimparam.nblks;
elseif isfield(Params.stimparam,'repetitions') 
    nblks2 = Params.stimparam.repetitions;
else
    nblks2=1;
end

if nblks < nblks2 && stimleninblk>1
    bincomplete =1;
    nstimlastblk = mod(nblks*nstiminblk,nstiminblk);
    nblks =ceil(nblks);
else
    bincomplete =0;
end

for iblk=1:nblks
    if iblk == nblks && bincomplete,
        stimleninblk = nstimlastblk;
    end

    for istim=1:stimleninblk
        stpoint = marker((iblk-1)*nstiminblk+istim)+1;          
        if bincomplete,
            edpoint = stpoint + 1;
        else                
            edpoint = stpoint + Params.stimparam.stim_samplesinNI-1;

        end

        cond = stimtime(stpoint);
        inxes = stpoint : edpoint;
        stimtime(inxes) = cond;

    end
end    


tlen=max(length(frame_start),length(stimtime));
timeinfo.frame_start=zeros(tlen,1);
timeinfo.stimtime=zeros(tlen,1);
timeinfo.frame_start(1:length(frame_start))=frame_start(:);
timeinfo.stimtime(1:length(stimtime))=stimtime(:);



msperframe = Params.msperframe;
timescale = msperframe/1000;
inx=find(timeinfo.frame_start);
a=timeinfo.stimtime(inx);
inxstart=find(diff(a)>0);
inxend=find(diff(a)<0);
if length(inxend)<length(inxstart)
    inxend(end+1)=length(a);
end
stimlen = length(inxstart)
stimtimes=cell(1,stimlen);
stimframes=cell(1,stimlen);
for ie=1:stimlen    
    stimtimes{ie}=[inxstart(ie) inxend(ie)]*timescale;
    stimframes{ie}=[inxstart(ie) inxend(ie)];
end
stimtimeinframe.stimtimes = stimtimes;
stimtimeinframe.stimframes = stimframes;
stimtimeinframe. stimlen = stimlen;
timeinfo.stimtimeinframe = stimtimeinframe;