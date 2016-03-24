%% cellbody identification

% gch_mc2d = reshape(gch_mc,[yd*xd td]);
% meanimg=reshape(mean(gch_mc2d,2),[yd xd]);
fn = 'D:\data_2photon\061713 M39 Ryan\TSeries-06172013-1334-021\AVG_MC_TSeries-06172013-1334-021.tif'
fn = 'D:\data_2photon\061713 M39 Ryan\TSeries-06172013-1334-022\AVG_MC_TSeries-06172013-1334-022.tif'
fn = 'D:\data_2photon\061713 M39 Ryan\TSeries-06172013-1334-023\AVG_MC_TSeries-06172013-1334-023.tif'

% fn='D:\data_2photon\061213 Tg1 nowhisker Ryan\TSeries-06122013-1541-016\AVG_MC_TSeries-06122013-1541-016.tif'
meanimg = double(imread(fn));
% figure; imagesc(meanimg); colormap('gray'); axis equal
para.CaInd ='GCamp6';
para.Hb=1; para.Lb=4;
para.sc=1;
% automated selection
[out]=autosel(meanimg, para);



template = out.template;
cci = out.cci;
cci_org = out.cciorg;
inxboundary = out.inxboundary;
inx_inroiwithboundary_org = out.inx_inroiwithboundary_org;
inxboundary_org = out.inxboundary_org;


% manual selection
inx=cell2mat(inx_inroiwithboundary_org);
meanimg1 = (meanimg- min(meanimg(:)))/(max(meanimg(:))-min(meanimg(:)));
fnsave = '0612-1541-016.mat';
baseimg = repmat(meanimg1,[1 1 3]);
Para.meanimg = meanimg;
Para.baseimg = baseimg;
Para.fnsave = fnsave
Para.figname = 'Manual Selection'
Para.mode = 'GCamp6'
select_cluster2(cci,Para)

%verification
% a=load(fnsave);
% cci = a.newcci;
figure;
img01=(meanimg-min(meanimg(:)))/(max(meanimg(:))-min(meanimg(:)));
img01 = repmat(img01,[1 1 3]);
hImg = overlayImage(img01, cci, [],'Discrete')
Nroi = length(unique(cci(:)))-1;
listroi = setdiff(unique(cci(:))',0);
for iroi= 1:Nroi
    inxs = find(cci(:)==listroi(iroi));
    [yc,xc] = ind2sub(size(cci),inxs(round(end/2)));
    text(xc,yc,num2str((iroi)),'Color','w');
end

