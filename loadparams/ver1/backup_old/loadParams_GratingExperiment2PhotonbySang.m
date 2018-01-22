function Params = loadParams_GratingExperiment2PhotonbySang(Params)

mainpath = Params.files.mainpath;
edrfn = Params.files.edrfn;
subpath_xml = Params.files.subpath_xml;
xml_fn = Params.files.xml_fn;
subpath_stim = Params.files.subpath_stim;
stim_fn1 = Params.files.stim_fn1;
ScreenRefreshrate = Params.ScreenRefreshrate;
if isfield(Params.files,'stim_fn2')
    stim_fn2 = Params.files.stim_fn2;
else
    stim_fn2 = '';
end


fullfn_edr = fullfile(mainpath, edrfn);
[edr_data, ~] = import_edr(fullfn_edr);
samplingint = diff(edr_data(:,1));
Params.samplingfreq_NI = 1/samplingint(1);

tinfo = (1:size(edr_data,1))/Params.samplingfreq_NI;
hfig = figure; plot(tinfo,edr_data(:,3));
title('mirro signal for identification of time period search')
xlabel('time (sec)'); ylabel('edr ch3')
uiwait(hfig);
prompt = {'Select the sample range tp search.'};
dlg_title = 'Timestamp for the first and the last samples';
num_lines = 1;
def = {['[0 ' num2str(tinfo(end)) ']']};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(answer)
    edrsamplerange= floor(str2num(answer{1})*Params.samplingfreq_NI)+1;    
    if edrsamplerange(2)>size(edr_data,1)
        edrsamplerange(2) = size(edr_data,1);
    end
    edr_data = edr_data(edrsamplerange(1):edrsamplerange(2),:);
    edr_data(:,1) = edr_data(:,1) - edr_data(1,1);
end

% figure; plot(edr_data(:,3))
fullfn_xml = fullfile(mainpath, subpath_xml, xml_fn);
xml_data = xml_parseany(fileread(fullfn_xml));

Params.max_frame = max(size(xml_data.Sequence{1,1}.Frame));
Params.Nframes = Params.max_frame;

lastframetime = str2double(xml_data.Sequence{1,1}.Frame{1,Params.max_frame}.ATTRIBUTE.relativeTime);
firstframetime = str2double(xml_data.Sequence{1,1}.Frame{1,1}.ATTRIBUTE.relativeTime);
Params.msperframe = 1000*( lastframetime - firstframetime )/(Params.max_frame-1);


Params.timeNI = edr_data(:,1);
mirror1 = edr_data(:,3);
mirror2 = edr_data(:,4);
Params.mirror1 = mirror1;
Params.mirror2 = mirror2;

Nsamples_NI = size(edr_data,1);
%     for ii=2:7
%     figure(ii); plot(edr_data(1:end,1), edr_data(1:end,ii))
%     end
%     legend('2','3','4','5','6','7');
%     ylim([-100 500])
% zmirror1 = mirror1 - mean(mirror1);
% tinx= find(edr_data(:,1)>12 & edr_data(:,1)<12.7);

% 
zmirror1= ztrans(mirror1,1);
zmirror2 = ztrans(mirror2,1);
dM = abs(diff(zmirror1))+abs(diff(zmirror2));
% 

bphoto = questdlg('Do you want to identify the timestamp of frames from photo diode?', ...
    'Frame timestamp indetification method', ...
    'Yes', 'No', 'Yes');

% 
% A= fftshift(fft(zmirror1));
% figure; plot(abs(A))
% figure; plot(mirror1-mean(mirror1))
if ~isempty(bphoto) && strcmp(bphoto,'Yes')
if strcmp(Params.scan_mode,'spiral')
% In my observation,when the mirror2 reached a lowest value,
% both mirrors settle on the gradual movements in the first frame.
% Therefore, I refer the frame start time point as the gradual start after the abrupt change of mirror2
% position.    
    [~, initguess_start] = max(dM(1:round(Nsamples_NI/4)));    
%     tinx =( initguess_start-round(0.01*Params.samplingfreq_NI)):(initguess_start+round(0.3*Params.samplingfreq_NI));
%     figure; plot(edr_data(tinx,1), [dM(tinx) zmirror1(tinx) zmirror2(tinx)],'.-')
%     figure; plot( [dM(tinx) zmirror1(tinx) zmirror2(tinx)],'.-')
%     figure; plot(zmirror2(initguess_start:initguess_start+5))

    clear mirror_frame_start mirror_start_time
    [~, offset] = min(zmirror2(initguess_start:initguess_start+5));
    ik = initguess_start + offset-1;
    mirror_frame_start(ik) = 10;
    mirror_start_time(1) = ik;
    
    estsamplesbtwframe = round(Params.msperframe/1000*Params.samplingfreq_NI);
    x = 2;
    ixs=[];
    while ik<Nsamples_NI
        
        istart = ik + round(estsamplesbtwframe*0.9);
        iend  = ik + round(estsamplesbtwframe*1.1);
        if iend> Nsamples_NI
            break;
        end
        ixs=[ixs;[istart iend]];

%         sdv = diff(dM(istart:iend));        
%         asdv = abs(sdv);
%         [~,inx1]=max(sdv);        
        
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
    
    if length(mirror_start_time)>= Params.Nframes
        ddf=diff(diff(mirror_start_time))*1000/Params.samplingfreq_NI;
        figure; plot(ddf); title('Difference of frame duration (msec)');
        bTimestampOfframe = questdlg('Do you accept the timestamp of frames identification?', ...
                         'Frame timestamp indetification', ...
                         'Yes', 'No', 'Yes');
        if strcmp(bTimestampOfframe,'Yes')             
            mirror_start_time = mirror_start_time(1:Params.Nframes);
            Params.mirror_start_time = mirror_start_time;  
            clear mirror_frame_start;
            mirror_frame_start (mirror_start_time)=10;
            Params.mirror_frame_start = mirror_frame_start;      
            Params.diff_mirror_start_time = diff(mirror_start_time);


            image_frame_duration = (mirror_start_time(size(mirror_start_time,2))-mirror_start_time(1))/(size(mirror_start_time,2)-1);
            for i = 1: length(mirror_start_time)
                image_frame_time(i) = round(mirror_start_time(i)+image_frame_duration/2);
                image_frame(image_frame_time(i)) = 1;
            end
            Params.frame_times = image_frame_time;
            Params.image_frame = image_frame;
        end
        

    elseif length(mirror_start_time)< Params.Nframes
        error('frame identification failed; Identified frames are less than the Nframes');
    end
else
    error('not implemented yet');
end

else
    [~, initguess_start] = max(dM(1:round(Nsamples_NI/4)));  
    bTimestampOfframe='No';
end
%% load scanning parameters
Params.date_time = xml_data.ATTRIBUTE.date;
numpar = length(xml_data.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key);
for ipa=1:numpar
    keyname = xml_data.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,ipa}.ATTRIBUTE.key;
    keyvalue = xml_data.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,ipa}.ATTRIBUTE.value;
    Params.extrascanparam.(keyname) = keyvalue;
end
%------- when photodiode resolution is not good enough, frame
%identification can be done via reading xml file and check whether xml file
%is correct.
if strcmp(bTimestampOfframe,'No')
    hfig=figure; 
    plot([zmirror1(1:initguess_start+1000) zmirror2(1:initguess_start+1000)]);
    legend('mirror-1','mirror-2');
    title('Re-identification of timestamp of the first frame');
    uiwait(hfig);
    prompt = {'Select the timestamp for the first frame start.'};
    dlg_title = 'Timestamp for the first frame';
    num_lines = 1;
    def = {num2str(initguess_start)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if ~isempty(answer)
        firstframeinx= str2num(answer{1});
        mirror_frame_start = zeros(size(mirror1));
        mirror_start_time = zeros(1,Params.Nframes);
        for iframe = 1 : Params.Nframes 
            ik = firstframeinx + round(str2num(xml_data.Sequence{1}.Frame{iframe}.ATTRIBUTE.relativeTime) * Params.samplingfreq_NI);
            mirror_frame_start(ik) = 10;    
            mirror_start_time(iframe) = ik;
        end
        
        
        ik2= ik-1000;
        if ik2<=0,
            ik2=1;
        end
        hfig= figure; 
        plot([mirror_frame_start(ik2:end) zmirror2(ik2:end)]);
        legend('From xml file','mirror2');
        title('Timestamp identification of frames');
        uiwait(hfig);
        bTimestampOfframe = questdlg('Do you accept the timestamp of frames identification?', ...
                         'Frame timestamp indetification', ...
                         'Yes', 'No', 'Yes');
                     
        if strcmp(bTimestampOfframe,'Yes') 
            mirror_start_time = mirror_start_time(1:Params.Nframes);
            Params.mirror_start_time = mirror_start_time;  
            clear mirror_frame_start;
            mirror_frame_start (mirror_start_time)=10;
            Params.mirror_frame_start = mirror_frame_start;      
            Params.diff_mirror_start_time = diff(mirror_start_time);


            image_frame_duration = (mirror_start_time(size(mirror_start_time,2))-mirror_start_time(1))/(size(mirror_start_time,2)-1);
            for i = 1: length(mirror_start_time)
                image_frame_time(i) = round(mirror_start_time(i)+image_frame_duration/2);
                image_frame(image_frame_time(i)) = 1;
            end
            Params.frame_times = image_frame_time;
            Params.image_frame = image_frame;
        else
            error('Identification of the timestamp for frames failure');
        end


    end
end

%% loading stimulus parameters and stimulus time from photo diode.
fullfn_stim1 = fullfile(mainpath, subpath_stim,stim_fn1);
load(fullfn_stim1)
paramconstants = stim.params.constants;
paramtrial = stim.params.trials;

stimparam.ScreenRefreshrate = ScreenRefreshrate;
stimparam.stim_samplesinNI = paramtrial(end).stimFrames*Params.samplingfreq_NI/ScreenRefreshrate;
stimparam.blank_samplesinNI = paramtrial(end).blankFrames*Params.samplingfreq_NI/ScreenRefreshrate;
stimparam.total_samplesinNI = paramtrial(end).stimulusTime*Params.samplingfreq_NI;
stimparam.repetitions = round(stimparam.total_samplesinNI/(stimparam.stim_samplesinNI + stimparam.blank_samplesinNI));
stimparam.contrast = paramtrial(end).contrast;
stimparam.orientation = paramtrial(end).orientation;
stimparam.spatialFreq = paramtrial(end).spatialFreq;
stimparam.tempoFreq = paramtrial(end).tempoFreq;



fullfn_stim2 = fullfile(mainpath, subpath_stim,stim_fn2);
if exist(fullfn_stim2)==2,
    Folder= fullfile(mainpath,subpath_stim);
    s = nex2mat(Folder,stim_fn2);
    timestamps = 1:stimparam.repetitions;
    s = setTimeStamp(s,timestamps);
    %retrieve the full set of stimulus events lookup table from s.
    %retrieve marker from nex file
    marker = s.nexData.markers{1};
    %stim variable set
    encode_vars = cell(1,length(marker.values));
    for j = 1 : length(marker.values)
        encode_vars{j} = marker.values{j}.name;
    end
    %keywords for event-sorting. full set of LUT returned when all variables
    %are true, i.e, '>0'
    event = struct;
    for j = 1 : length(encode_vars)
        event(j).type = encode_vars{j};
        event(j).string = '>0';
        event(j).operator = '&';
    end
    %retrive the full set of stim-event-lookup table.
    [t_SETS,StimEventLUT] = sortStimEvent(s,event);
    cond.t_SETS = t_SETS;
    cond.StimEventLUT =StimEventLUT;
    cond.event = event;
end
stimparam.cond = cond;




% figure; plot(edr_data(:,5))
% load the photo-diode signal
PDsig = edr_data(:,5); 

smoothPDsig = smooth (PDsig, 2*Params.samplingfreq_NI);
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
% %get the stimulus onet times from teh photodiode trace
% sorted_start = sort(flatPDsig(round(interval(1)*Params.samplingfreq_NI):round(interval(2)*Params.samplingfreq_NI)),'descend');
% % start_thr = mean(sorted_start(1:round(0.1*length(sorted_start))));
% thresh_PD=sorted_start(1)*2
% tstamp =round( interval(2)*Params.samplingfreq_NI)+1;

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
tstamp0=tstamp;
tstamp = tstamp -100;
if tstamp <0,
    tstamp =1;    
end
initstep = tstamp0-tstamp-10;

clear stim_onset_tstampinNI onset_tstamp_vector blank_onset_tstampinNI

kk = 1;
jj = 1;
for iblk = 1:stimparam.repetitions   
    if iblk==1,
        tstamp = tstamp + initstep;
    else
        tstamp = tstamp + round(stimparam.blank_samplesinNI*0.9);
    end
    if tstamp<length(flatPDsig)
        while  flatPDsig(tstamp)<thresh_PD 
            if tstamp<length(flatPDsig)
                tstamp=tstamp+1;
            end
        end 
        stim_onset_tstampinNI(1,kk)=tstamp;
        onset_tstamp_vector(1,tstamp) = 10;
        kk = kk + 1;  

        tstamp = tstamp + round(stimparam.stim_samplesinNI*0.9);
        while  flatPDsig(tstamp)<thresh_PD 
            if tstamp<length(flatPDsig)
                tstamp=tstamp+1;
            end
        end 
        blank_onset_tstampinNI(1,jj)=tstamp;
        onset_tstamp_vector(1,tstamp) = -10;
        jj = jj + 1;        
    end
end

len=length(onset_tstamp_vector);
len2 = length(flatPDsig)
onset_tstamp_vector2 =[onset_tstamp_vector'; zeros(len2-len,1)];
hfig=figure; plot(Params.timeNI(1:len2),[onset_tstamp_vector2*100 flatPDsig(1:len2)])

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