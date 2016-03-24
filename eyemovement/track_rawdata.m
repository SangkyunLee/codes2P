function [out, MSK]= track_rawdata(MOV,opts, InitPar)
% function track_rawdata(MOV,opts)
% opts.winsize=10;
% opts.thr =30
% opts.bdisp =false;


winsize = opts.winsize; 
thr = opts.thr;

%% select mask
if nargin<3
    SEL.thr = thr ;
    SEL.op ='>';
    InitPar = gen_manual_pupilmask(MOV,1, SEL);
end

if nargout>2
    MSK = uint8(zeros(size(MOV,1)*size(MOV,2),size(MOV,3)));    
end
%% create clusters
Nframe = size(MOV,3);
CM = cell(Nframe,1);
CM{1} = InitPar.CM(1:2);
PIX = cell(Nframe,1);
PIX{1} = {InitPar.pix};
parfor iimg = 1:Nframe
    if iimg==1,
        continue;
    end   
    
    MSK1 = zeros(size(MOV,1),size(MOV,2));
    fmap = zeros(size(MOV,1),size(MOV,2));
    fmap(MOV(:,:,iimg)<thr)=1;

    [Li, num] = bwlabel(double(fmap),8);
    for ic=1:num
        cpix = find(Li(:)==ic);
        if length(cpix)<100
            Li(cpix)=0;
            MSK1(cpix)=128;
        end
    end
    CID = setdiff(unique(Li(:)),0);
    Nc = length(CID);        
    CM0=zeros(2,Nc);
    PIX0=cell(1,Nc);
    for ic = 1 : Nc
        Lx = zeros(size(Li));
        Cpix1 = find(Li(:)==CID(ic));
        Lx(Cpix1) = ic;
        [X0,Y0]=cmass(Lx);
        CM0(:,ic)=[X0 Y0];            
        PIX0{ic} = single(Cpix1);
    end

    PIX{iimg}=PIX0;
    CM{iimg}=CM0;
    MSK(:,iimg)=uint8(MSK1(:));
end
pause(1);

%% find the closest cluster as the pupil  
sPIX = cell(Nframe,1);
CInx = zeros(Nframe,1);
Cinfo = zeros(2,Nframe);
Cinfo(:,1)=CM{1};
dists = zeros(Nframe,1);
sPIX{1} = PIX{1}{1};
if opts.bdisp
    figure; 
end



for iimg = 2 : Nframe    
    CM1 =CM{iimg};
    Nc=size(CM1,2);
    dist=zeros(winsize,Nc);
    for iwin = 1: winsize
        if iimg<=iwin,
            continue;
        end
        dist(iwin,:) = sqrt(sum(bsxfun(@minus,CM1(1:2,:),Cinfo(1:2,iimg-iwin)).^2,1));
    end
    dist = mean(dist,1);
    
    A = zeros(size(MOV,1),size(MOV,2));
    [~, inx]=min(dist);
    dists(iimg)=dist(inx);
    CInx(iimg)=inx;
    Cinfo(:,iimg)=CM{iimg}(:,inx);
    sPIX{iimg} = PIX{iimg}{inx};
    for ic = 1: length(dist)    
        if ic == inx
            A(PIX{iimg}{inx})=255;
            MSK(PIX{iimg}{ic},iimg) = 255;
        else
            MSK(PIX{iimg}{ic},iimg)=128;
        end
    end
        
    
    if opts.bdisp
        imagesc(A);
        title(num2str(iimg));
        pause(0.002);
    end
end

out.sPIX = sPIX;
out.CM = Cinfo;
MSK = reshape(MSK, size(MOV));

