function Params = loadParams_Imagepresentation2(Params)

mainpath = Params.files.mainpath;
edrfn = Params.files.edrfn;
subpath_xml = Params.files.subpath_xml;
xml_fn = Params.files.xml_fn;
subpath_stim = Params.files.subpath_stim;
stim_fn1 = Params.files.stim_fn1;
ScreenRefreshrate = Params.ScreenRefreshrate;


fullfn_edr = fullfile(mainpath, edrfn);
[edr_data, h] = import_edr(fullfn_edr);
samplingint = diff(edr_data(:,1));
Params.samplingfreq_NI = 1/samplingint(1);


fullfn_xml = fullfile(mainpath, subpath_xml, xml_fn)
xml_data = xml_parseany(fileread(fullfn_xml));

Params.max_frame = max(size(xml_data.Sequence{1,1}.Frame));
Params.Nframes = Params.max_frame;

lastframetime = str2num(xml_data.Sequence{1,1}.Frame{1,Params.max_frame}.ATTRIBUTE.relativeTime);
firstframetime = str2num(xml_data.Sequence{1,1}.Frame{1,1}.ATTRIBUTE.relativeTime);
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
zmirror1 = mirror1 - mean(mirror1);
% tinx= find(edr_data(:,1)>12 & edr_data(:,1)<12.7);

% 
zmirror1= ztrans(mirror1,1);
zmirror2 = ztrans(mirror2,1);
dM = abs(diff(zmirror1))+abs(diff(zmirror2));
% 


% 
% A= fftshift(fft(zmirror1));
% figure; plot(abs(A))
% figure; plot(mirror1-mean(mirror1))

if strcmp(Params.scan_mode,'spiral')
% In my observation,when the mirror2 reached a lowest value,
% both mirrors settle on the gradual movements in the first frame.
% Therefore, I refer the frame start time point as the gradual start after the abrupt change of mirror2
% position.    
    [~, initguess_start] = max(dM(1:round(Nsamples_NI/4)));    
    tinx =( initguess_start-round(0.01*Params.samplingfreq_NI)):(initguess_start+round(0.3*Params.samplingfreq_NI));
%     figure; plot(edr_data(tinx,1), [dM(tinx) zmirror1(tinx) zmirror2(tinx)],'.-')
%     figure; plot( [dM(tinx) zmirror1(tinx) zmirror2(tinx)],'.-')
%     figure; plot(zmirror2(initguess_start:initguess_start+5))

    clear mirror_frame_start mirror_start_time
    [~, offset] = min(zmirror2(initguess_start:initguess_start+5))
    ik = initguess_start + offset-1;
    mirror_frame_start(ik) = 10;
    mirror_start_time(1) = ik;
    
    estsamplesbtwframe = round(Params.msperframe/1000*Params.samplingfreq_NI);
    x = 2;
    ixs=[];
    while ik<Nsamples_NI
        
        istart = ik + round(estsamplesbtwframe*0.6);
        iend  = ik + round(estsamplesbtwframe*1.4);
        if iend> Nsamples_NI
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
    
    if length(mirror_start_time)>= Params.Nframes
        ddf=diff(diff(mirror_start_time));
        figure; plot(ddf);
        
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
        

    elseif length(mirror_start_time)< Params.Nframes
        error('frame identification failed; Identified frames are less than the Nframes');
    end
    
    
end

%% load scanning parameters
Params.date_time = xml_data.ATTRIBUTE.date;
numpar = length(xml_data.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key);
for ipa=1:numpar
    keyname = xml_data.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,ipa}.ATTRIBUTE.key;
    keyvalue = xml_data.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,ipa}.ATTRIBUTE.value;
    Params.extrascanparam.(keyname) = keyvalue;
end

%% loading stimulus parameters and stimulus time from photo diode.
fullfn_stim1 = fullfile(mainpath, subpath_stim,stim_fn1);
load(fullfn_stim1)
paramconstants = stim.params.constants;
paramtrial = stim.params.trials;


stimparam.ScreenRefreshrate = ScreenRefreshrate;
stimparam.stim1_samplesinNI = (1/paramconstants.stimframerate) * Params.samplingfreq_NI;
stimparam.stimblk_samplesinNI = paramconstants.stimtimeinblk * Params.samplingfreq_NI;
chksumforstimtype = stimparam.stimblk_samplesinNI-stimparam.stim1_samplesinNI;

stimparam.blankblk_samplesinNI = paramconstants.blanktimeinblk * Params.samplingfreq_NI;
stimparam.nblks = paramconstants.nblks;
stimparam.nstiminblk = round(stimparam.stimblk_samplesinNI/stimparam.stim1_samplesinNI);
stimparam.imgOrder = stim.params.trials(end).imgOrder;
stimparam.Imgmatfile = stim.params.constants.Imgmatfile;


% figure; plot(edr_data(:,5))
% load the photo-diode signal
PDsig = edr_data(:,5); 

smoothPDsig = smooth (PDsig, 2*Params.samplingfreq_NI);
flatPDsig = PDsig -mean(PDsig);

len=round(length(flatPDsig)*0.2);
h=figure; plot(edr_data(1:len,1),flatPDsig(1:len)); xlabel('time (sec)');
waitfor(h)

prompt = {'Time interval'};
dlg_title = 'Select time interval for the inital threshold for the photo-diode';
num_lines = 1
def = {'[10 11]'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
interval = eval(answer{1});

%get the stimulus onet times from teh photodiode trace
sorted_start = sort(flatPDsig(round(interval(1)*Params.samplingfreq_NI):round(interval(2)*Params.samplingfreq_NI)),'descend');
% start_thr = mean(sorted_start(1:round(0.1*length(sorted_start))));
thresh_PD=sorted_start(1)*2

tstamp =round( interval(2)*Params.samplingfreq_NI)+1;
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
nstiminblk = stimparam.nstiminblk;
if nstiminblk ==1,
    for iblk = 1:stimparam.nblks   
        if iblk==1,
            tstamp = tstamp + initstep;
        else
            tstamp = tstamp + round(stimparam.blankblk_samplesinNI*0.9);
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

            tstamp = tstamp + round(stimparam.stimblk_samplesinNI*0.9);
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
else
    nstiminblk2 = stimparam.nstiminblk+1;
    for iblk = 1:stimparam.nblks        
        for is = 1 : nstiminblk2
            if iblk==1,
                tstamp = tstamp + initstep;
            else
                tstamp = tstamp + round(stimparam.stim1_samplesinNI*0.9) ;
            end
            
            if tstamp<length(flatPDsig)
                while  flatPDsig(tstamp)<thresh_PD 
                    if tstamp<length(flatPDsig)
                        tstamp=tstamp+1;
                    end
                end 

                if is == nstiminblk2
                    blank_onset_tstampinNI(1,jj)=tstamp;
                    onset_tstamp_vector(1,tstamp) = -10;
                    jj = jj + 1;
                else
                    stim_onset_tstampinNI(1,kk)=tstamp;
                    onset_tstamp_vector(1,tstamp) = 10;
                    kk = kk + 1;  
                end            
            end

        end
        tstamp =  tstamp + round(stimparam.blankblk_samplesinNI*0.9) ;
    end
end

len=length(onset_tstamp_vector);
len2 = length(flatPDsig)
onset_tstamp_vector2 =[onset_tstamp_vector'; zeros(len2-len,1)];
hfig=figure; plot(edr_data(1:len2),[flatPDsig(1:len2) onset_tstamp_vector2(1:len2)*20])
onset_tstamp_vector = onset_tstamp_vector2;

%%

Nstiminblk = stimparam.nstiminblk;
Nimgpresented = length(stim_onset_tstampinNI)
Nblks = ceil(Nimgpresented/Nstiminblk);
a = zeros(size(Params.timeNI));

inxblank = 0;
for iimg=1:Nimgpresented
    inxstart=stim_onset_tstampinNI(iimg);
    if mod(iimg,Nstiminblk)==0,
        inxblank = inxblank + 1; 
        inxend=blank_onset_tstampinNI(inxblank)-1;
    else
        iimg2= iimg+1;
        if iimg2>length(stim_onset_tstampinNI)
            break;
        end
        inxend=stim_onset_tstampinNI(iimg2)-1;
    end    
    a(inxstart:inxend)=stimparam.imgOrder(iimg);
end



%%




stimparam.stim_onset_tstampinNI = stim_onset_tstampinNI;
stimparam.blank_onset_tstampinNI = blank_onset_tstampinNI;
stimparam.onset_tstamp_vector = onset_tstamp_vector;

Params.stim_onset_tstampinNI = stim_onset_tstampinNI;
Params.blank_onset_tstampinNI = blank_onset_tstampinNI;
Params.onset_tstamp_vector = onset_tstamp_vector;
Params.ScreenRefreshrate = ScreenRefreshrate;
Params.stimparam = stimparam;
Params.stim_tstamp_vector = a;

uiwait(hfig);