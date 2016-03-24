%% cellbody identification

% gch_mc2d = reshape(gch_mc,[yd*xd td]);
% meanimg=reshape(mean(gch_mc2d,2),[yd xd]);
fn = 'D:\data_2photon\061213 Tg1 nowhisker Ryan\accum12-17\Mean12-17.tif'
meanimg = double(imread(fn));
% figure; imagesc(meanimg); colormap('gray'); axis equal
para.CaInd ='GCamp6';
para.Hb=1; para.Lb=3; 
[out]=autosel(meanimg, para);



template = out.template;
template_org = out.templateorg;
cci = out.cci;
cci_org = out.cciorg;
inxboundary = out.inxboundary;
inx_inroiwithboundary_org = out.inx_inroiwithboundary_org;
inxboundary_org = out.inxboundary_org;


% selection ROI 
inx=cell2mat(inx_inroiwithboundary_org);
meanimg1 = (meanimg- min(meanimg(:)))/(max(meanimg(:))-min(meanimg(:)));
fnsave = '';
baseimg = repmat(meanimg1,[1 1 3]);
select_cluster2(meanimg,cci,baseimg,fnsave,'Manual Removal')


% manual seperation
% close all
cellids=[14 23 26 49 60 71 77 112];
% manualselection
[inx_inroiwithboundary,inxboundary]=manualselection(cellids, template_org, inx_inroiwithboundary_org,inxboundary_org);

% visual confirmation of the manual separation
CCInew=zeros(size(template_org));
for cellid=1:length(inxboundary)
    inxin=setdiff(inx_inroiwithboundary{cellid},inxboundary{cellid});
    CCInew(inxin)=cellid;
end
figure; imagesc(CCInew)
hFig=figure;
select_cluster(hFig,template_org,CCInew,cell2mat(inxboundary))