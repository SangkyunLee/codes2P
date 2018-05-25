function Params = loadParams_onlineRetinotopy(Params)

[stimparam_file,PathName] = uigetfile('*.edr','Select the edr file');
Params.files.mainpath=PathName;
Params.files.edrfn = stimparam_file;


mainpath = Params.files.mainpath;
edrfn = Params.files.edrfn;

% path_stim = Params.files.path_stim;
stim_fn1 = Params.files.stim_fn1;
ScreenRefreshrate = Params.ScreenRefreshrate;

fullfn_edr = fullfile(mainpath, edrfn);
[edr_data, h] = import_edr(fullfn_edr);
samplingint = diff(edr_data(:,1));
Params.samplingfreq_NI = 1/samplingint(1);


Params.timeNI = edr_data(:,1);
mirror1 = edr_data(:,3);
mirror2 = edr_data(:,4);
Params.mirror1 = mirror1;
Params.mirror2 = mirror2;

Nsamples_NI = size(edr_data,1);
zmirror1 = mirror1 - mean(mirror1);


% 
zmirror1= ztrans(mirror1,1);
zmirror2 = ztrans(mirror2,1);
dM = abs(diff(zmirror1))+abs(diff(zmirror2));
% 


if strcmp(Params.scan_mode,'spiral')
% In my observation,when the mirror2 reached a lowest value,
% both mirrors settle on the gradual movements in the first frame.
% Therefore, I refer the frame start time point as the gradual start after the abrupt change of mirror2
% position.    
    [~, initguess_start] = max(dM(1:round(Nsamples_NI/4)));    
    tinx =( initguess_start-round(0.01*Params.samplingfreq_NI)):(initguess_start+round(0.3*Params.samplingfreq_NI));

    clear mirror_frame_start mirror_start_time
    [~, offset] = min(zmirror2(initguess_start:initguess_start+5))
    ik = initguess_start + offset-1;
    mirror_frame_start(ik) = 10;
    mirror_start_time(1) = ik;
    
    estsamplesbtwframe = round(Params.msperframe/1000*Params.samplingfreq_NI);
    x = 2;
    ixs=[];
    while ik<Nsamples_NI-1
        
        istart = ik + round(estsamplesbtwframe*0.6);
        iend  = ik + round(estsamplesbtwframe*1.4);
        if iend> Nsamples_NI-1
            iend  = Nsamples_NI-1
        end
        if iend<=istart, 
            break; 
        end
        ixs=[ixs;[istart iend]];

        sdv = diff(dM(istart:iend));        
        asdv = abs(sdv);
        [~,inx1]=max(sdv);        
        
        amp =abs(zmirror2(istart:iend));
%         figure; plot(amp,'.-')
        A = [amp(1:end-5) amp(2:end-4) amp(3:end-3) amp(4:end-2) amp(5:end-1) amp(6:end)];
        dstdA = diff(std(A,0,2));
%         figure; plot(dstdA)
        [~, inx2] = max(dstdA);
        b = zmirror2(istart:iend);    
        if inx2+8>length(b)
            [~, inx3] = min(b(inx2:length(b)));
        else
            [~, inx3] = min(b(inx2:inx2+8));
        end
        ik = istart + (inx2-1)+(inx3-1);
        mirror_frame_start(ik) = 10;
        mirror_start_time(x) = ik;
        x = x + 1;                      
    end
    Nframe_identified = length(mirror_start_time);
    
    if Nframe_identified>= Params.Nframes
        ddf=diff(diff(mirror_start_time));
%         figure; plot(ddf);
        
        Params.mirror_start_time = mirror_start_time(1:Params.Nframes);
        clear mirror_frame_start;
        Params.mirror_frame_start(mirror_start_time)=10;        
        Params.diff_mirror_start_time = diff(mirror_start_time);
        
        
        image_frame_duration = (mirror_start_time(size(mirror_start_time,2))-mirror_start_time(1))/(size(mirror_start_time,2)-1);
        for i = 1: length(mirror_start_time)
            image_frame_time(i) = round(mirror_start_time(i)+image_frame_duration/2);
            image_frame(image_frame_time(i)) = 1;
        end
        Params.frame_times = image_frame_time;
        Params.image_frame = image_frame;
        
        
    elseif Nframe_identified < Params.Nframes
        msg = sprintf('Frame identification failed; Identified frames(%d) are less than the Nframes(%d).',Nframe_identified, Params.Nframes);
        warndlg(msg);
        Params.Nframes = length(mirror_start_time);
    end
    
    
end



%% loading stimulus parameters and stimulus time from photo diode.
fullfn_stim1 = fullfile(stim_fn1);
load(fullfn_stim1)
paramconstants = stim.params.constants;
paramtrial = stim.params.trials;

stimparam.ScreenRefreshrate = ScreenRefreshrate;
stimparam.stim_samplesinNI = paramtrial(end).stimFrames*Params.samplingfreq_NI/ScreenRefreshrate;
stimparam.blank_samplesinNI = paramtrial(end).blankFrames*Params.samplingfreq_NI/ScreenRefreshrate;
stimparam.total_samplesinNI = paramtrial(end).stimulusTime*Params.samplingfreq_NI;
stimparam.repetitions = round(stimparam.total_samplesinNI/(stimparam.stim_samplesinNI + stimparam.blank_samplesinNI));




% figure; plot(edr_data(:,5))
% load the photo-diode signal
PDsig = edr_data(:,5); 

flatPDsig = PDsig -mean(PDsig);

len=round(length(flatPDsig)*0.2);
h=figure; plot(edr_data(1:len,1),flatPDsig(1:len)); xlabel('time (sec)')
waitfor(h)

% prompt = {'Time interval'};
% dlg_title = 'Select time interval for the inital threshold for the photo-diode';
% num_lines = 1
% def = {'[10 11]'};
% answer = inputdlg(prompt,dlg_title,num_lines,def);
% interval = eval(answer{1});
% 
% %get the stimulus onet times from teh photodiode trace
% sorted_start = sort(flatPDsig(round(interval(1)*Params.samplingfreq_NI):round(interval(2)*Params.samplingfreq_NI)),'descend');
% % start_thr = mean(sorted_start(1:round(0.1*length(sorted_start))));
% thresh_PD=sorted_start(1)*2

% i =round( interval(2)*Params.samplingfreq_NI)+1;

prompt = {'Threshold','Time start'};
dlg_title = 'Input a starting point for the photo-diode';
num_lines = 2
def = {'300','5'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
thresh_PD = eval(answer{1});
inxt = eval(answer{2});

tstamp =round( inxt*Params.samplingfreq_NI)+1;

while flatPDsig(tstamp)<thresh_PD
    tstamp=tstamp+1;
end

start = tstamp;
clear stim_onset_tstampinNI onset_tstamp_vector blank_onset_tstampinNI
stim_onset_tstampinNI(1,1) = start;
onset_tstamp_vector(1,start) = 10;


j = start;
for i = 1:stimparam.repetitions    
    if stimparam.blank_samplesinNI>0
        j = j + round(stimparam.stim_samplesinNI*0.9) ;
        if j<length(flatPDsig)
            while flatPDsig(j)<thresh_PD
                if j<length(flatPDsig)
                    j=j+1;
                end
            end        
        end
        blank_onset_tstampinNI(1,i)=j;
        onset_tstamp_vector(1,j) = -10;
    end
    if i<stimparam.repetitions
        if stimparam.blank_samplesinNI>0
            j = j + round(stimparam.blank_samplesinNI*0.9) ;
        else
            j = j + round(stimparam.stim_samplesinNI*0.9) ;
        end
        if j<length(flatPDsig)
            while flatPDsig(j)<thresh_PD
                if j<length(flatPDsig)
                    j=j+1;
                end
            end
        end
        stim_onset_tstampinNI(1,i+1)=j;
        onset_tstamp_vector(1,j) = 10;
    end    
end

len=length(onset_tstamp_vector);
len2 = length(flatPDsig)
onset_tstamp_vector2 =[onset_tstamp_vector'; zeros(len2-len,1)];
hfig=figure; plot(edr_data(1:len2),[onset_tstamp_vector2*100 flatPDsig(1:len2)])

stimparam.stim_onset_tstampinNI = stim_onset_tstampinNI;

if exist('blank_onset_tstampinNI')
    stimparam.stim_duration_tstampinNI = blank_onset_tstampinNI - stim_onset_tstampinNI;
    stimparam.blank_onset_tstampinNI = blank_onset_tstampinNI;
else
    dff = stim_onset_tstampinNI(2:end) - stim_onset_tstampinNI(1:end-1);
    stimparam.stim_duration_tstampinNI = [dff round(mean(dff(2:end)))];
end
stimparam.onset_tstamp_vector = onset_tstamp_vector2;

Params.stim_onset_tstampinNI = stim_onset_tstampinNI;
Params.stim_duration_tstampinNI = stimparam.stim_duration_tstampinNI;
if exist('blank_onset_tstampinNI')    
    Params.blank_onset_tstampinNI = blank_onset_tstampinNI;
else
end
Params.onset_tstamp_vector = onset_tstamp_vector2;
Params.ScreenRefreshrate = ScreenRefreshrate;
Params.stimparam = stimparam;

uiwait(hfig);
