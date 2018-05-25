clear all
close all
Nframe=3500
samplingrate = 7.18 %hz


Params.scan_mode = 'spiral'
Params.Nframes = Nframe;
Params.msperframe = 1000/samplingrate;


stimparam.CanvasSizeDeg = 50;
stimparam.BarSizeDeg = 10;
stimparam.BarMovStepDeg = 3;
stimparam.BarFramesLoc = 1*60;
stimparam.blankFramesLoc = 2*60;
stimparam.ScreenRefreshrate =60;
stimparam.dir = [0 90]
stimparam.nrepetition =5;
stimparam.blankFrames = 2*stimparam.ScreenRefreshrate;
Params = loadParams_onlineBarRetinotopy(Params, stimparam);


fext='tif';
deli='_';
Nframe=Params.Nframes;
data = load_data(fext,deli,Nframe);   

% fext='zip';
% [fndata1,path_data1] = uigetfile(['*.' fext],'Select the image file');
% data = loadTseries([path_data1 filesep fndata1]);


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
data = data(:,1:Nframe);
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


plot_barRetinotopy_results(dFF,Params)
