clear all
close all

%---------- setting loadin parameters in lab1
% Nframe=2900
% samplingrate = 8.9%hz
% 
% Params.files.stim_fn1='Retinotopy_4son_X5_Y3.mat'
% Params.scan_mode = 'spiral'
% Params.ScreenRefreshrate =60; % 60Hz\
% Params.Nframes = Nframe;
% Params.msperframe = 1000/samplingrate;


%---------- setting loadin parameters in NB208
Nframe=2000
samplingrate = 6.3%hz
Params.Channel.inx_mrr1 = 3;
Params.Channel.inx_mrr2 = 3;
Params.Channel.inx_photo = 2;

Params.files.stim_fn1='Retinotopy_4son_X5_Y3.mat'
Params.scan_mode = 'spiral_NB208'
Params.ScreenRefreshrate =60; % 60Hz\
Params.Nframes = Nframe;
Params.msperframe = 1000/samplingrate;


%---- stimulus condition specification

Nrepetition =5;
Cond.Ncond =15;
Cond.Nx =5;
Cond.Ntrials =Cond.Ncond*Nrepetition;

% Cond.Ny =3;

%----loading parameters 
Params = loadParams_onlineRetinotopy(Params)
pause(1);
%---data collection
fext='ome.tif';
deli='_';
Nframe=Params.Nframes;
data = load_data(fext,deli,Nframe);   

% fext='zip';
% [fndata1,path_data1] = uigetfile(['*.' fext],'Select the image file');
% data = loadTseries([path_data1 filesep fndata1]);

%--- ROI selection 
roigui(data)

%% mask mask and display
dsize=size(data);
meanimg = mean(data,3);
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

data = reshape(data,[dsize(1)*dsize(2) dsize(3)]);
Fi=[];
listroi = setdiff(unique(MASK(:))',0);
for iroi=listroi
    inxvxs=find(MASK(:)==iroi);    
    avgsig=mean(data(inxvxs,:),1);
    Fi=[Fi avgsig(:)];
end

period=round(20*1000/Params.msperframe);
Fb=[];
for frame_i = period:(Nframe-period)
    around_i = sort(Fi((frame_i-period+1):frame_i+period,:));        
    Fb(frame_i,:) = mean(around_i(1:round(length(around_i)*0.2),:));    
end

Fb(1:period,:) =repmat(Fb(period,:),[period 1]);
Fb(Nframe-period+1:Nframe,:) = repmat(Fb(Nframe-period,:),[period 1]);    
dFF=(Fi-Fb)./Fb;
plotopt=struct
plot_retinotopy_results(dFF,Params,Cond,[])

