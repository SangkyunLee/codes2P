function [out]=autosel(meanimg, para)
Hb=para.Hb; % difference of Gaussian filter
Lb=para.Lb; % difference of Gaussian filter

% calcium indicatior; OGB, GCamp6
if isfield(para,'CaInd')
    CaInd = para.CaInd; 
    if strcmp(CaInd,'OGB')
        if isfield(para,'cs')
            cs=para.cs;  % morphological filter           
        else
            cs =[];
        end
    elseif strcmp(CaInd,'GCamp6')
        
        cs = [];
    end
else
    CaInd = [];
end

if isfield(para,'thr')
    thr = para.thr; %threshold for mask 
else
    thr = 0;
end

if isfield(para,'minpix_thr_cluster')
    minpix_thr_cluster = para.minpix_thr_cluster; %threshold after clustering
else
    minpix_thr_cluster = 1;
end
if isfield(para,'sc')    
    sc=para.sc; % scale factor for image enlargement
else
    sc = 1;
end
if isfield(para,'fnsave')
    fnsave = para.fnsave;
else
    fnsave =  ['temp-' datestr(now, 'mmddyyyy-HHMMSS')];
end

meanimg = imresize(meanimg,sc);
%------ for spiral scan ----------
inx0 = find(meanimg(:)==0);
a=zeros(size(meanimg));
a(inx0) =1;
a1 = filterDoG(a, Hb, Lb);
inx0 = [inx0; find(a1(:)~=0)];
inx0 = unique(inx0);
%--------------------------------


template = filterDoG(meanimg, Hb, Lb);
template(inx0)=0;
inxn0 = setdiff([1:prod(size(template))],inx0);
template(inxn0) = (template(inxn0) -min(template(inxn0)))/(max(template(inxn0))-min(template(inxn0)));
template_3d = repmat(template,[1 1 3]);

ztem = zeros(size(template));
ztem(inxn0) = ztrans(template(inxn0),2);
ztem2d = reshape(ztem,size(template));
% figure; imagesc(ztem2d);

%---------------------------------
baseimg = template_3d;
fmap = ztem2d;
% fnsave='aa.mat'
autosel_GUI(meanimg,fmap,baseimg,'Continuous',fnsave,'1st Threshold in single voxels');
THR=load(fnsave);
if ~THR.bvalid
    error('no selection of the threhold values');
end
thr =THR.thr;
%--------------------------------

ztem = sign( ztem -thr);
ztem(find(ztem<0))=0;
mask =reshape(ztem,[size(template)]);
% figure; imagesc(mask); axis equal;

if strcmp(CaInd,'OGB') 
    if ~isempty(cs)
        se = strel('disk',cs);
        mask = imopen(mask,se);
    end
    
    [cci,roinum] = spm_bwlabel(mask,6);
    h=[0 1 0; 1 0 1; 0 1 0];
    inxboundary={};
    for iroi=1:roinum    
        inx=find(cci(:)==iroi);
        a=zeros(size(cci));
        a(inx)=1;
        Y=filter2(h,a);
        Y(find(Y(:)==4))=0;
        inxboundary{iroi} = find(Y(:)~=0)';
    end

    % save cellidentify meanimg template_org mask ztem
    inx=cell2mat(inxboundary);
    inxin=setdiff([1:prod(size(mask))],inx);
else
    [cci,roinum] = spm_bwlabel(mask,6);
    inxboundary = {};
end

% thresholding and rorder
newcci = cci;
newid =1;
for ic=1:roinum
   inxp_clust = find(cci==ic);
   if length(inxp_clust) < minpix_thr_cluster
       newcci(inxp_clust)=0;
   else
       newcci(inxp_clust) = newid;
       newid = newid + 1;
   end
end
cci = newcci;
%---------------------------------
meanimg = (meanimg- min(meanimg(:)))/(max(meanimg(:))-min(meanimg(:)));
autosel_GUI(meanimg,cci,repmat(meanimg,[1 1 3]), 'Discrete',fnsave,'2nd Thresh. in Clust.');
THR=load(fnsave);
if ~THR.bvalid
    error('no selection of the threhold values');
end
thr2 =THR.thr
%--------------------------------
minpix_thr_cluster =thr2;
newcci = cci;
newid =1;
for ic=1:roinum
   inxp_clust = find(cci==ic);
   if length(inxp_clust) < minpix_thr_cluster
       newcci(inxp_clust)=0;
   else
       newcci(inxp_clust) = newid;
       newid = newid + 1;
   end
end
cci = newcci;

%% scale back to the original domain

sc2=1/sc;
if sc2~=1
    mask_org=imresize(sign(cci),sc2,'nearest');
   [cciorg,roinum] = spm_bwlabel(mask_org,6);
else
    cciorg = cci;
    mask_org = mask;
    roinum = length(unique(cci))-1;
end
if strcmp(CaInd,'OGB') 
    h=[0 1 0; 1 0 1; 0 1 0];
    inxboundary_org={};
    inx_inroiwithboundary_org ={};
    % figure;
    for iroi=1:roinum    
        inx=find(cciorg(:)==iroi);
        a=zeros(size(mask_org));
        a(inx)=1;
        Y=filter2(h,a);
        Y(find(Y(:)>3))=0;
        inxboundary_org{iroi} = find(Y(:)~=0)';

        b=zeros(size(mask_org));
        b(find(Y(:)~=0)')=1;
        inx_inroiwithboundary_org{iroi}= find(sign(a+b)==1)';
    end
else
    inxboundary_org ={};
    for iroi=1:roinum                
        inx_inroiwithboundary_org{iroi}= find(cciorg==iroi)';
    end
end




out.template = template;
% out.templateorg = template_org;
out.cci = cci;
out.cciorg = cciorg;
out.inxboundary = inxboundary;

out.inxboundary_org = inxboundary_org;
out.inx_inroiwithboundary_org = inx_inroiwithboundary_org;


