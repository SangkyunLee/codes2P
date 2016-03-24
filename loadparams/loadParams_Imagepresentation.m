function Params = loadParams_ImagePresentation(Params)
% function Params = loadParams_ImagePresentation(Params)
% written by Sangkyun Lee, v1 2013-x-x
% Params.files.stim_fn1='Retinotopy_4son_X5_Y3.mat'
% Params.scan_mode = 'spiral' | 'spiral_NB208'
% Params.ScreenRefreshrate =60; % 60Hz\
% Params.Nframes = Nframe;
% Params.msperframe = 1000/samplingrate;
% 
% modified by Sangkyun Lee, v2 2015-08-10 for NB208 in the neurosensory
% building
% Params.Channel.inx_mrr1 = 3;
% Params.Channel.inx_mrr2 = 3;
% Params.Channel.inx_photo = 2;
% Params.files.stim_fn1='Retinotopy_4son_X5_Y3.mat'
% Params.scan_mode = 'spiral' | 'spiral_NB208'
% Params.ScreenRefreshrate =60; % 60Hz\
% Params.Nframes = Nframe;
% Params.msperframe = 1000/samplingrate;
% 


if isfield(Params,'Channel')
    inx_mirror1 = Params.Channel.inx_mrr1;
    inx_mirror2 = Params.Channel.inx_mrr2;
    inx_photodiode = Params.Channel.inx_photo;
else
    if strcmp(Params.scan_mode,'spiral')
        inx_mirror1=3;
        inx_mirror2=4;
        inx_photodiode=5;
    elseif strcmp(Params.scan_mode,'spiral_NB208')
        inx_mirror1=3;
        inx_mirror2=3;
        inx_photodiode=2;
    elseif strcmp(Params.scan_mode,'resonant_NB208')
        inx_mirror1=3;
        inx_mirror2=3;
        inx_photodiode=2;
    else
        error('not implemented yet');
    end
    
end
   
    

mainpath = Params.files.mainpath;
DAQfn = Params.files.DAQfn;

path_stim = Params.files.subpath_stim;
stim_fn1 = Params.files.stim_fn1;
ScreenRefreshrate = Params.ScreenRefreshrate;

fullfn_DAQ = fullfile(mainpath, DAQfn);

xml_fn = Params.files.xml_fn;
subpath_xml = Params.files.subpath_xml;
fullfn_xml = fullfile(mainpath, subpath_xml, xml_fn)
xml_data = xml_parseany(fileread(fullfn_xml));
Params.Nframes = max(size(xml_data.Sequence{1,1}.Frame));

lastframetime = str2num(xml_data.Sequence{1,1}.Frame{1,Params.Nframes}.ATTRIBUTE.relativeTime);
firstframetime = str2num(xml_data.Sequence{1,1}.Frame{1,1}.ATTRIBUTE.relativeTime);
Params.msperframe = 1000*( lastframetime - firstframetime )/(Params.Nframes-1);



if strcmp(DAQfn(end-3:end),'.EDR')

    [DAQ_data, h] = import_edr(fullfn_DAQ);
    samplingint = diff(DAQ_data(:,1));
    Params.samplingfreq_NI = 1/samplingint(1);


    Params.timeNI = DAQ_data(:,1);
    mirror1 = DAQ_data(:,inx_mirror1);
    mirror2 = DAQ_data(:,inx_mirror2);
    Params.mirror1 = mirror1;
    Params.mirror2 = mirror2;

    Nsamples_NI = size(DAQ_data,1);
    zmirror1 = mirror1 - mean(mirror1);
    % 
    zmirror1= ztrans(mirror1,1);
    zmirror2 = ztrans(mirror2,1);
    dM = abs(diff(zmirror1))+abs(diff(zmirror2));
elseif strcmp(DAQfn(end-3:end),'.csv') %-------- this code for NB208 in neurosensory building
    % This function is too slow
    %M = importdata(fullfn_DAQ);
    %DAQ_data = M.data;
    fid = fopen(fullfn_DAQ);
    C1 = textscan(fid, '%s %s%d, %s%d, %s%d, %s%d\n');
    DAQ_data = textscan(fid, '%f, %f, %f, %f, %f','CollectOutput', 1);
    fclose(fid)    
    DAQ_data = DAQ_data{1};    
    DAQ_data(:,1) = DAQ_data(:,1)/1000;
    
    Params.timeNI = DAQ_data(:,1);
    samplingint = diff(DAQ_data(:,1));
    Params.samplingfreq_NI = 1/samplingint(1);
    
    
    mirror1 = DAQ_data(:,inx_mirror1);
    mirror2 = DAQ_data(:,inx_mirror2);
    Params.mirror1 = mirror1;
    Params.mirror2 = mirror2;

    Nsamples_NI = size(DAQ_data,1);
    zmirror1 = mirror1 - mean(mirror1);
    % 
    zmirror1= ztrans(mirror1,1);
    zmirror2 = ztrans(mirror2,1);
    dM = abs(diff(zmirror1));
    dM = [dM(:); 0];
    
    
end

switch Params.scan_mode
    case {'spiral','spiral_NB208'}
        test=1
    otherwise
        test=0;
end
        

bphoto = questdlg('Do you want to identify the timestamp of frames from photo diode?', ...
    'Frame timestamp indetification method', ...
    'Yes', 'No', 'Yes');

% 
% A= fftshift(fft(zmirror1));
% figure; plot(abs(A))
% figure; plot(mirror1-mean(mirror1))
if ~isempty(bphoto) && strcmp(bphoto,'Yes')

    
    switch Params.scan_mode
        case {'spiral','spiral_NB208'}
            
           
            % In my observation,when the mirror2 reached a lowest value,
            % for NB208 in the neurosensory building, mirror1 is the same
            % as mirror2
            % both mirrors settle on the gradual movements in the first frame.
            % Therefore, I refer the frame start time point as the gradual start after the abrupt change of mirror2
            % position.    
            if strcmp(Params.scan_mode,'spiral_NB208')
                DAQsamplingperiod = Params.timeNI(2)-Params.timeNI(1);
                Nsample = round(5/DAQsamplingperiod); %search 5seconds of the beginnig of voltage acqusition
                [~, initguess_start] = max(dM(1:Nsample));    
            else                
                [~, initguess_start] = max(dM(1:round(Nsamples_NI/4)));    
            end
            % tinx =( initguess_start-round(0.01*Params.samplingfreq_NI)):(initguess_start+round(0.3*Params.samplingfreq_NI));
            % figure; plot(edr_data(tinx,1), [dM(tinx) zmirror1(tinx) zmirror2(tinx)],'.-')
            % figure; plot( [dM(tinx) zmirror1(tinx) zmirror2(tinx)],'.-')
            % figure; plot(zmirror2(initguess_start:initguess_start+5))

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
                    mirror_start_time = mirror_start_time(1:Params.Nframes+1);
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
            
        %------ resonant galvanometer. 
        % In the following code, it is required to consider the averaged
        % frames
        case {'resonant', 'resonant_NB208'}
            hfig=figure;    
            plot(Params.timeNI(1:5000),[dM(1:5000) zmirror1(1:5000)]);
            legend('Mirror Diff','mirror-1');
            title('Iidentification of timestamp of the first frame');
            uiwait(hfig);
            prompt = {'Select the threshold for frame start.'};
            dlg_title = 'Threshold';
            num_lines = 1;
            def = {'1'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            Thr = str2num(answer{1});
            guessTnewfrm = find(dM>Thr);
            initguess_start = guessTnewfrm(1);

        
            ik = initguess_start;
            mirror_frame_start = zeros(size(zmirror1));
            clear mirror_start_time;
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
                [~,inx]=max(dM(istart:iend));
                
                ik = istart + inx-1;
                if dM(ik)>Thr
                    mirror_frame_start(ik) = 10;
                    mirror_start_time(x) = ik;
                    x = x + 1;  
                end
                  
            end

            if length(mirror_start_time)>= Params.Nframes
                ddf=diff(diff(mirror_start_time))*1000/Params.samplingfreq_NI;
                figure; plot(ddf); title('Difference of frame duration (msec)');
                bTimestampOfframe = questdlg('Do you accept the timestamp of frames identification?', ...
                                 'Frame timestamp indetification', ...
                                 'Yes', 'No', 'Yes');
                if strcmp(bTimestampOfframe,'Yes')             
                    mirror_start_time = mirror_start_time(1:Params.Nframes+1);
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
                bphoton='No'
            end
                

            
    
        otherwise
            error('not implemented yet');
    end

end

if strcmp(bphoto,'No')  
    if strcmp(Params.scan_mode,'spiral_NB208')
        DAQsamplingperiod = Params.timeNI(2)-Params.timeNI(1);
        Nsample = round(2/DAQsamplingperiod); %search 5seconds of the beginnig of voltage acqusition
        [~, initguess_start] = max(dM(1:Nsample));    
    elseif strcmp(Params.scan_mode,'spiral')                
        [~, initguess_start] = max(dM(1:round(Nsamples_NI/4)));    
    else
        intiguess_start=0;
    end 
    bTimestampOfframe='No';
    
    fullfn_xml = fullfile(Params.files.mainpath,Params.files.subpath_xml,Params.files.xml_fn);
    xml_data = xml_parseany(fileread(fullfn_xml));
    Params.date_time = xml_data.ATTRIBUTE.date;
%     numpar = length(xml_data.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key);
%     for ipa=1:numpar
%         keyname = xml_data.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,ipa}.ATTRIBUTE.key;
%         keyvalue = xml_data.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,ipa}.ATTRIBUTE.value;
%         Params.extrascanparam.(keyname) = keyvalue;
%     end
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
end



%% loading stimulus parameters and stimulus time from photo diode.
fullfn_stim1 = fullfile(mainpath,path_stim,stim_fn1);
load(fullfn_stim1)
paramconstants = stim.params.constants;
paramtrial = stim.params.trials;

stimparam.ScreenRefreshrate = ScreenRefreshrate;
stimparam.stim1_samplesinNI = (1/paramconstants.stimframerate) * Params.samplingfreq_NI;
stimparam.stimblk_samplesinNI = paramconstants.stimtimeinblk * Params.samplingfreq_NI;
stimparam.blankblk_samplesinNI = paramconstants.blanktimeinblk * Params.samplingfreq_NI;
stimparam.nblks = paramconstants.nblks;
stimparam.nstiminblk = round(stimparam.stimblk_samplesinNI/stimparam.stim1_samplesinNI);
stimparam.imgOrder = stim.params.trials(end).imgOrder;
stimparam.Imgmatfile = stim.params.constants.Imgmatfile;


% figure; plot(edr_data(:,5))
% load the photo-diode signal
PDsig = DAQ_data(:,inx_photodiode); 

flatPDsig = PDsig -mean(PDsig);
flatPDsig = abs(PDsig);

len=round(length(flatPDsig)*0.2);
h=figure; plot(DAQ_data(1:len,1),flatPDsig(1:len)); xlabel('time (sec)')
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
onset_tstamp_vector = zeros(1, length(zmirror1));
stim_onset_tstampinNI(1,1) = start;
onset_tstamp_vector(1,start) = 10;


j = start;
for i = 1 : stimparam.nblks
    if stimparam.blankblk_samplesinNI>0
        j = j + round(stimparam.stimblk_samplesinNI*0.9) ;
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
    if i<stimparam.nblks
        if stimparam.blankblk_samplesinNI>0
            j = j + round(stimparam.blankblk_samplesinNI*0.9) ;
        else
            j = j + round(stimparam.stimblk_samplesinNI*0.9) ;
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
hfig=figure; plot(DAQ_data(1:len2),[onset_tstamp_vector2/5*(max(flatPDsig(1:len2))) flatPDsig(1:len2)])

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

