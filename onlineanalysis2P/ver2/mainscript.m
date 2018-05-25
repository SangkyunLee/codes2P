clear all
close all

%%------- retiontopy 
Params.files.mainpath='Z:\Sang\020518_thy1_0819_1\2pimage\TSeries-02052018-0737-005';

Nframe=2300;
samplingrate = 7.2034; %hz
Params.Channel.inx_frtrg = 0;
Params.Channel.inx_photo = 1;
Params.Channel.nch =4;
Params.scan_mode = 'spiral';
Params.Nframes = Nframe;
Params.msperframe = 1000/samplingrate;

Params.files.stim_mainpath=[];
Params.files.stim_fn1='\\slcompss\slee\codes2P\onlineanalysis2P\ver2\Retinotopy_5X3Y_2on2off_5rep.mat';
Nrep =5;
Cond.Ncond =15;
Cond.Nx=5;
Cond.condseq = repmat(1:Cond.Ncond,[1 Nrep]);


%% ------------- orientation tuning
% 
% Params.files.mainpath='D:\packages\onlineanalysis2P\test_data\TSeries-12262017-1444-010';
% Nframe=2000;
% samplingrate = 7.2034; %hz
% Params.Channel.inx_frtrg = 1;
% Params.Channel.inx_photo = 2;
% Params.Channel.nch =2;
% Params.scan_mode = 'spiral';
% Params.Nframes = Nframe;
% Params.msperframe = 1000/samplingrate;
% 
% Params.files.stim_mainpath=[];
% Params.files.stim_fn1='D:\packages\onlineanalysis2P\ver2\Grating_12dir_p5secon_1p5secoff.mat';
% Cond.Ncond =12;
% Cond.condseq = [];
% 
% 
% %% ------------- size tuning
% 
% Params.files.mainpath='D:\packages\onlineanalysis2P\test_data\TSeries-01222018-0853-001';
% Nframe=1300;
% samplingrate = 7.2034; %hz
% Params.Channel.inx_frtrg = 1;
% Params.Channel.inx_photo = 0;
% Params.Channel.nch =2;
% Params.scan_mode = 'spiral';
% Params.Nframes = Nframe;
% Params.msperframe = 1000/samplingrate;
% 
% Params.files.stim_mainpath=[];
% Params.files.stim_fn1='D:\packages\onlineanalysis2P\test_data\CrossOrientationExperiment.mat';
% Cond.Ncond =12;
% Cond.condseq = [];


%%

%----loading parameters 
Params = loadParams(Params);
if strcmp(Params.VStimparam.expType,'GratingExperiment2PhotonbySang')
    Cond.condseq = Params.VStimparam.DIOValue;
end
timeinfo= gen_VStimtime(Params, Cond);


%---data collection
% loading sequential data
fext='tif';%'ome.tif';
deli='_';
Nframe=Params.Nframes;
tic
Tseries = load_data(fext,deli,Nframe);   
c=toc
%loading stack data
% fext='tif';
% [fndata1,path_data1] = uigetfile(['*.' fext],'Select the image file');
% Tseries = loadTseries([path_data1 filesep fndata1]);
% Tseries = Tseries(:,:,1:Nframe);

%--- ROI selection 
roigui(Tseries)

%% mask mask and display
dsize=size(Tseries);
meanimg = mean(Tseries,3);
meanimg1 = (meanimg- min(meanimg(:)))/(max(meanimg(:))-min(meanimg(:)));
baseimg = repmat(meanimg1,[1 1 3]);
MASK = zeros(dsize(1), dsize(2));
for iroi=1:length(ROI_list)
    inxs=ROI_list(iroi).pixel_list;
    MASK(inxs)=iroi;
end
figure; 
subplot(121);
imagesc(meanimg1); colormap('gray'); axis equal; xlim([1 size(meanimg1,2)]); ylim([1 size(meanimg1,1)])
subplot(122);
overlayImage(baseimg,MASK, [],'Discrete');
haxis = gca;
display_ROInum(haxis,MASK);


%--- extracting signals from individual ROIs and plot the results according to the given conditions

Tseries = reshape(Tseries,[dsize(1)*dsize(2) dsize(3)]);

listroi = setdiff(unique(MASK(:))',0);
Fi=zeros(dsize(3), length(listroi));
for i = 1: length(listroi)
    roi = listroi(i);
    inxvxs=find(MASK(:)==roi);    
    avgsig=mean(Tseries(inxvxs,:),1);
    Fi(:,i)= avgsig(:);
end

%%-- Fb calculation with temporal smoothing time window size = bgtime
bgtime =5;
period=round(bgtime*1000/Params.msperframe); % time period for background activity
Fb=[];
for frame_i = period:(Nframe-period)
    around_i = sort(Fi((frame_i-period+1):frame_i+period,:));        
    Fb(frame_i,:) = mean(around_i(1:round(length(around_i)*0.2),:));    
end

Fb(1:period,:) =repmat(Fb(period,:),[period 1]);
Fb(Nframe-period+1:Nframe,:) = repmat(Fb(Nframe-period,:),[period 1]);    
dFF=(Fi-Fb)./Fb;



%--------------
plotRet; % plot online retinotopy results
% plot_ORItuning; % plot orientation tuning results
