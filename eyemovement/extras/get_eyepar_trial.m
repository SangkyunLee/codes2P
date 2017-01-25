function [XY, R, dP, PP, vfinx] = get_eyepar_trial(Params,tinfo)
% function [XY, R, dP, PP, vfinx] = get_eyepar_trial(Params,tinfo)
% This function return the pupil center and radius sorted by trials
% XY, R, dP: cell(nevt,1), nevt: no trials
% dP: dif. of t and t+1
% PP: pupil parameters from all frames recorded
% vfinx: videoframe index sorted by stimulus trial
% Params: loaded by load_data and get_eyeparfn(T,scan,mp_os)
% tinfo: timeinfo loaded by gen_stimtime_xxx(Params);
% 08-23-2016, Sangkyun Lee

%- These predefined parameters can be changed.
prestim = 0.1; % prestimtime for eyetracking
poststim = 0.1; % post-stimtime for eyetracking

%-------------------------------------------

F = Params.files;
if ~isfield(F,'eyeparfn')
    error('Params.files does not contain eyeparfn');
else
    eyeparfn = F.eyeparfn;
end
mp = F.mainpath;
eyeparpath = fullfile(mp,'matlab','eyepar');
eyepar = load(fullfile(eyeparpath,eyeparfn));

fit = eyepar.fitInfo; 
info = eyepar.searchInfo;
Xoff = info.Xr(1);
Yoff = info.Yr(1);
nf = length(fit);
PP=zeros(nf,3);
for i = 1 : nf
    p = fit(i).par;
    if ~isempty(p)
        PP(i,:) = p;
    end
end
dPP = PP(2:end,:)-PP(1:end-1,:);
dPP = [zeros(1,3); dPP];



eye_dir = fullfile(mp,F.eyemovdir);
rn = str2double(F.eyemovfn); 
[daqdata,vidtime,~] = load_eyedaq(eye_dir,rn);
[stimstart_eyedaq, ~] = identify_stimstart(daqdata(:,1), daqdata(:,2),2 ,1);


% align eyemovietime and stim time from the frist stim start
sf  = Params.samplingfreq_NI;


stimtime = tinfo.stimtime;
dst = diff(stimtime);
stim_tstamp = [find(dst>0) find(dst<0)]/sf;
stim_rt = stim_tstamp - stim_tstamp(1); % relative stim time from stim start
eye_rt = vidtime-stimstart_eyedaq;  % relative eyemovie time from stim start



nevt = size(stim_rt,1);
inval = zeros(nevt,1);
XY = cell(nevt,1);
R = cell(nevt,1);
dP = cell(nevt,1);
vfinx = cell(nevt,1); % video frameinx in absolute time
for it = 1 : nevt
    st = stim_rt(it,1) - prestim;
    et = stim_rt(it,2) + poststim;
    
    inx = (eye_rt>st & eye_rt<et);
    PP0 = PP(inx,:);

    if any(PP0(:,3)==0)
        inval(it)=1;
    else
        XY{it} = [PP0(:,1)+Xoff-1 PP0(:,2)+Yoff-1];
        R{it} = PP0(:,3);
        dP{it} = dPP(inx,:);
        vfinx{it} = find(inx);
   
    end
    
end
