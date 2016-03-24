function [Y others]=data_sort(data,spec,bquiet, Params)
% function [Y others]=data_sort(data,spec)
%
% INPUT:
% frames = spec.frames;
% dataType = spec.dataType;
% nCell = spec.nCell;
% head motion in timeseries
% spec.motion.tmotion; 
% Threshold for trial rejection for big motion
% spec.motion.motionthr; 
% or 
% spec.motion.tmotion1(motionthr1) & spec.motion.tmotion2(motionthr2)
%
% OUTPUT:
% Y (Trial x frames x cells)
% others.stiminx
% others.stimlen
% others.events

% tmotion1 & tmmotion2 added 6/16/2014
% modified by sangkyun lee 04/15/2015 for separating argument for param

if nargin<3
    bquiet=false;
end
if nargin<4
    for idata=1:length(data)
        Params(idata) = data(idata).Params;
        data(idata).Params = struct([]);
    end
end

    
frames = spec.frames;
dataType = spec.dataType;
nCell = spec.nCell;
if isfield(spec,'motion') & isstruct(spec.motion)
    if isfield(spec.motion,'tmotion')
        tmotion = spec.motion.tmotion;
        motionthr = spec.motion.motionthr;
        Bigmotion = tmotion>motionthr; 
    elseif isfield(spec.motion,'tmotion1') && isfield(spec.motion,'motionthr1')
        tmotion1 = spec.motion.tmotion1;
        motionthr1 = spec.motion.motionthr1;
        Bigmotion = tmotion1 > motionthr1; 
        if isfield(spec.motion,'tmotion2') && isfield(spec.motion,'motionthr2')
            tmotion2 = spec.motion.tmotion2;
            motionthr2 = spec.motion.motionthr2;
            Bigmotion = eval(['(tmotion1 >motionthr1)' spec.motion.op ' (tmotion2 >motionthr2)']); 
        end
    else
        error('no motion parameter specified');
    end
    if isfield(spec.motion, 'bdisp') && spec.motion.bdisp
        figure; 
        plot([tmotion1; tmotion2; Bigmotion]');
        title('motion');
        legend('motion in pix','diff','Bigmotion');
    end
    
    bmotion = true;
else
    bmotion = false;
end


stimlen=zeros(1,length(data));
stiminx = cell(1,length(data));
for idata=1:length(data);   
    
    
    if isfield(data,'timeinfo') && ~isempty(data.timeinfo)
        timeinfo = data.timeinfo;
    elseif isfield(Params,'timeinfo')  && ~isempty(Params.timeinfo)
        timeinfo = Params.timeinfo;
    else
        if  strcmp(Params(idata).files.stim_fn1,'GratingExperiment2Photon.mat')
            timeinfo = gen_stimtime_GratingExperiment2Photon(Params(idata));
        elseif strcmp(Params(idata).files.stim_fn1,'GratingExperiment2PhotonbySang.mat')
            timeinfo = gen_stimtime_GratingExperiment2PhotonbySang(Params(idata));
        else
            error('not coded')
        end
    end
    stimtime = timeinfo.stimtime;
    frame_start = timeinfo.frame_start;
    inxf = find(frame_start==10);     

    if length(inxf)<Params(idata).Nframes 
        if ~bquiet
            msg = sprintf('Identified frames(%d) is less than the recorded frames(%d)', length(inxf),Params(idata).Nframes);
            warndlg(msg);
        end
        Nframes = length(inxf);
    else
        Nframes = Params(idata).Nframes;
    end
    inxf = inxf(1:Nframes);  
    inxin = length(stimtime)>=inxf;
    stimsel = stimtime(inxf(inxin));                
    stiminx{idata} = stimsel;

    stimlen(idata)=length(setdiff(unique(stimsel(:)),0));
end
others.stiminx = stiminx;
others.stimlen =stimlen;


%% resort data with stim onset
Y=cell(1,length(data));
events =cell(1,length(data));
for idata=1:length(data)
   
    sig = data(idata).(dataType);
    if nCell==size(sig,1)
        sig = sig';
    elseif nCell == size(sig,2)
    else
        error('data dimension mismatch');
    end
    Nframes = Params(idata).Nframes;
        

    stimonset = find(diff(stiminx{idata})~=0)+1;        
    trialdur = diff(stimonset);
    stimimginx = stiminx{idata}(stimonset);
    inxn0 = find(stimimginx~=0);
    stimonset = stimonset(inxn0);
    trialdur = trialdur(inxn0(1:end-1));
    trialdur(end+1) = round(mean(trialdur));
    stimimginx = stimimginx(inxn0); 
    
    % take out some trials of which windows lie outside of collected
    % frames.
    cut_pre=zeros(1,length(stimonset));
    cut_post=zeros(1,length(stimonset));
    if bmotion
        cut_motion = zeros(1,length(stimonset));
    end
    
    for itrial=1:length(stimonset)
        if any((stimonset(itrial)+frames)<1),
            cut_pre(itrial)=1;
        end
        if any((stimonset(itrial)+frames)>Nframes),
            cut_post(itrial)=1;
        end      
        if bmotion
            inxframe = stimonset(itrial)+[-5:trialdur(itrial)+5];
            if inxframe(1)>0 & inxframe(end)<=Nframes
                cut_motion(itrial) = any(Bigmotion(inxframe));
            end
        end
    end
    
    if bmotion
        keepinx = ~(cut_pre | cut_post | cut_motion);    
    else
        keepinx = ~(cut_pre | cut_post);    
    end
    stimonset = stimonset(keepinx);
    stimimginx = stimimginx(keepinx);
    
    events{idata} = stimimginx;
    sY=zeros(length(stimimginx),length(frames),nCell);
    for iwin=1:length(frames)        
        a = sig(stimonset+frames(iwin),:);
        a=reshape(a,[size(a,1) 1 size(a,2)]);
        sY(:,iwin,:)=a; 
    end
    Y{idata}=sY;
    seltrial{idata} = find(keepinx);
end
others.events = events(:); 
others.seltrial = seltrial;

