%% identify cell body from a scan
clear all

% first scan
mainpath='D:\data_2photon\130827BL\RF_0827BL_S1_131007\data\'
subpath='TSeries-10072013-1549-001'
fn='MC_TSeries-10072013-1549-001_Cycle00001_CurrentSettings_Ch2.zip';
finfo=fullfile(mainpath,subpath,fn);
Img =loadTseries(finfo);
Img = single(Img);
img1 = mean(Img,3);
% roigui(Img)
% save('ROI_1549-001.mat','ROI_list');
roigui(Img,ROI_list);
%% Check the ROIs from all scans

iscan=2
subpath=sprintf('TSeries-10072013-1549-%03d',iscan);
fn=sprintf('AVG_MC_TSeries-10072013-1549-%03d_Cycle00001_CurrentSettings_Ch2.tif',iscan);
finfo=fullfile(mainpath,subpath,fn);
Img =loadTseries(finfo);
Img = single(Img);
img1 = mean(Img,3);
M=ROIoverlap(img1,ROI_list);
% roigui(Img,ROI_list)
% save('ROI_1549-001.mat','ROI_list');

%% ROI sort by spatial location
nROI= length(ROI_list);
[m n] = size(img1);
Pos = zeros(2, nROI);
for iroi = 1:nROI
    Pos(:,iroi) = ROI_list(iroi).centerPos(:);
end
inx = sub2ind([m,n],round(Pos(1,:)),round(Pos(2,:)));
[mv newinx]=sort(inx,'ascend');

clear newROI;
for iroi=1:nROI    
    newROI(iroi) = ROI_list(newinx(iroi));
    newROI(iroi).name = num2str(iroi);
end
roigui(Img,newROI);

%%
close all
clear newROIs
iscan=17
subpath=sprintf('TSeries-10072013-1549-%03d',iscan);
fn=sprintf('AVG_MC_TSeries-10072013-1549-%03d_Cycle00001_CurrentSettings_Ch2.tif',iscan);
finfo=fullfile(mainpath,subpath,fn);
Img =loadTseries(finfo);
Img = single(Img);
img1 = mean(Img,3);
M=ROIoverlap(img1,newROI);
uiwait(M.fig);

fnsave = sprintf('ROI_1549-%03d.mat',iscan)
ROI = newROIs;
save(fnsave,'ROI');
% a=zeros(200);
% for i=1:88
% a(ROI(i).neuropilarea)=i;
% end

