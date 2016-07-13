function [mask M ] = est_localmotion_parameters(refimg,meanimg, images, fn_alignmask, winsize, bgpu, bglobal)
% [mask M ] = est_localmotion_parameters(refimg,meanimg, images, fn_alignmask, winsize)
% winsize: Estimation of motion parameter after temporal averaging or
% weighted average when winsize is a vector, winsize should be always odd
% number or a odd size of vector
% bgpu: this option is added for GPU processing. 2015-07-21, Sangkyun Lee


if nargin<5,
    fn_alignmask=[];
    winsize=1;
    bgpu = false;
    bglobal = false;
end

if nargin<6,
    bgpu = false;
    bglobal=false;
end
if nargin <7
    bglobal = false;
end


data= images;


if ~bglobal
    
    inxn0 = find(meanimg(:)~=0);
    % mask1 = zeros(size(meanimg));
    % mask1(inxn0)=meanimg(inxn0);
    mask2 = zeros(size(meanimg));
    masktmp= var(data,0,3);
    mask2(inxn0) = masktmp(inxn0);

    meanimg1 = (meanimg- min(meanimg(:)))/(max(meanimg(:))-min(meanimg(:)));
    fmap = mask2;
    baseimg = repmat(meanimg1,[1 1 3]);
    fnsave ='test.mat'

    autosel_GUI(meanimg1,fmap,baseimg,'Continuous',fnsave,'1st Threshold in single voxels');
    Thr=load(fnsave);
    Thr= Thr.thr;
    mask2(mask2<Thr)=0;
    figure; imagesc(log(1+mask2)); colormap('gray')

    meanimg1 = (meanimg- min(meanimg(:)))/(max(meanimg(:))-min(meanimg(:)));
    fnsave = 'temp.mat';
    baseimg = repmat(meanimg1,[1 1 3]);
    Para.meanimg = meanimg;
    Para.baseimg = baseimg;
    Para.fnsave = fnsave;
    Para.figname = 'Manual Selection';
    Para.mode = 'GCamp6';
    select_regions(mask2,Para);

    %% select ROIs
    maskroi=load(fnsave);
    hfig=figure; 
    subplot(121);
    imagesc(meanimg1); colormap('gray'); axis equal; xlim([1 size(meanimg1,2)]); ylim([1 size(meanimg1,1)])
    subplot(122);
    overlayImage(baseimg,maskroi.newcci, [],'Discrete');
    haxis = gca;
    display_ROInum(haxis,maskroi.newcci);
    %         uiwait(hfig);
    %         display('identified')

%     %% individual motion correction in selected areas 
%     if  length(winsize)>1
%         winsize=winsize(:);
%         hwin = floor(length(winsize)/2);
%         winsize = winsize/sum(winsize);
%         win3d= reshape(winsize, [1 1 length(winsize)]);
%         win3d = repmat(win3d, [size(data,1) size(data,2)]);
%     else    
%         hwin = floor(winsize/2);
%     end
%     cci=maskroi.newcci;
%     mask=zeros(size(cci));
% 
%     T = size(data,3);
% 
%     f_BP = filterDoG(refimg,3,15);        
%     listROIs = maskroi.selectedROIs;
%     nroi = length(listROIs);

else
    maskroi.newcci = ones(size(data,1),size(data,2));
    maskroi.selectedROIs = 1;
    
end


if  length(winsize)>1
    winsize=winsize(:);
    hwin = floor(length(winsize)/2);
    winsize = winsize/sum(winsize);
    win3d= reshape(winsize, [1 1 length(winsize)]);
    win3d = repmat(win3d, [size(data,1) size(data,2)]);
else    
    hwin = floor(winsize/2);
end
cci=maskroi.newcci;
mask=zeros(size(cci));

T = size(data,3);

f_BP = filterDoG(refimg,3,15);        
listROIs = maskroi.selectedROIs;
nroi = length(listROIs);


%------------------
outs_1=zeros(nroi,4,T);


for iimg=1:T
    if length(winsize)>1
        iimgs = (iimg-hwin:iimg+hwin);
        iimgs = iimgs(iimgs>0 & iimgs<=T);
        tmpImg = data(:,:,iimgs);
        if size(tmpImg,3)<length(winsize)
            g = mean(data(:,:,iimgs),3);
        else
            g = sum(tmpImg.*win3d,3);
        end
    else
        iimgs = (iimg-hwin:iimg+hwin);
        iimgs = iimgs(iimgs>0 & iimgs<=T);
        g = mean(data(:,:,iimgs),3);
    end
    
    g_BP = filterDoG(g,3,15);
    outinfo=[];
    for inxroi=listROIs
        inxrg = find(maskroi.newcci==inxroi);
        [ys, xs] = ind2sub(size(maskroi.newcci),inxrg);
        yr = min(ys): max(ys); xr = min(xs): max(xs);
        if bgpu
            %tic
            fsub = gpuArray(f_BP(yr,xr));
            gsub = gpuArray(g_BP(yr,xr));
            f1 = (fft2(fsub));            
            g1 = (fft2(gsub));           
            
            [output] = dftregistration(f1,g1,3, bgpu);
            output = gather(output);
            %toc;
        
        else
            %tic 
            fsub=f_BP(yr,xr);
            gsub=g_BP(yr,xr);   
            [output] = dftregistration(fft2(fsub),fft2(gsub),3,0);
            %toc
        end
        outinfo = [outinfo; output];        
    end
    outs_1(:,:,iimg) = outinfo;    
end



err = reshape((outs_1(:,1,:)),[nroi T]);
hfig=figure;
for iroi = 1:nroi
    subplot(nroi,1,iroi);
    plot(err(iroi,:));    
    title([num2str(iroi) ', err: ' num2str(round(100*mean(err(iroi,:)))) ]);
end

yshift = reshape((outs_1(:,3,:)),[nroi T]);
hy=figure;
for iroi=1:nroi
    subplot(nroi,1,iroi);
    plot(yshift(iroi,:));
    title(['Yshfit: ' num2str(iroi)]);
%         ylim([-10 10]);
end

xshift = reshape((outs_1(:,4,:)),[nroi T]);
hx=figure;
for iroi=1:nroi
    subplot(nroi,1,iroi);
    plot(xshift(iroi,:));
    title(['Xshfit: ' num2str(iroi)])
%             ylim([-5 5])
end
% uiwait(hfig);


hmsg= msgbox(['mean error:' num2str(mean(err,2)')]);
uiwait(hmsg);
if exist('hfig','var') 
    try
        close(hfig);
    catch err        
        disp(err.message);
    end
end
if exist('hx','var')   
    try
        close(hx);
    catch err        
        disp(err.message);
    end
end
if exist('hy','var')   
    try
        close(hy);
    catch err        
        disp(err.message);
    end    
end

if nroi>1
    prompt = {'Select ROI numbers'};
    dlg_title = 'ROI numbers for motion correction';
    num_lines = 1;
    def = {'[1]'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    a=answer{1};
    eval(['selroi=' a]);
else
    selroi=1;
end


ymargin=5; xmargin=5;
for iroi =selroi
    inxc = find(cci(:)==iroi);
    [ys xs]= ind2sub(size(cci),inxc);
    cy=mean(ys); cx=mean(xs);
    radiusy = (max(ys)-cy);
    radiusx = (max(xs)-cx);
    scaley = (1: ymargin);
    scaley = scaley/length(scaley);        
    scalex = (1: xmargin);
    scalex = scalex/length(scalex);
    
    
    ygrid = cy-(radiusy+ymargin) : cy+(radiusy+ymargin);
    xgrid = cx-(radiusx+xmargin) : cx+(radiusx+xmargin);
    ylen = length(ygrid);
    xlen = length(xgrid);
    
    ay = [scaley(:); ones(ylen-2*length(scaley),1); flipud(scaley(:))];
    ax = [scalex(:); ones(xlen-2*length(scalex),1); flipud(scalex(:))];
    a = ay*ax';
    inxy = find(ygrid>0 & ygrid<size(cci,1));
    inxx = find(xgrid>0 & xgrid<size(cci,2));
    
%     a=hamming(radiusy*2+1)*hamming(radiusx*2+1)';
    
    mask(ygrid(inxy),xgrid(inxx))=a(inxy,inxx);
end
figure; imagesc(mask);

if length(selroi)==1,
    M=squeeze(outs_1(selroi,:,:));
else
    M=[];
end
if ~isempty(fn_alignmask)
save(fn_alignmask,'mask','M');
end

