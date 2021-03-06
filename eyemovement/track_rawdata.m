function [out, MSK]= track_rawdata(MOV,opts, InitPar, newCM,htwin)
% function track_rawdata(MOV,opts)
% opts.winsize=10;
% opts.thr =30
% opts.bdisp =false;


winsize = opts.winsize; 
thr = opts.thr;
if nargin<5
    htwin=0;
end
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
    if (iimg-htwin)<1 || (iimg+htwin)>Nframe
        a=MOV(:,:,iimg);
        fmap(a(:)<thr)=1;
    else
        a=mean(MOV(:,:,iimg-htwin:iimg+htwin),3);
        fmap(a(:)<thr)=1;
    end

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
%     if iimg==717,
%         keyborad;
%     end
    CM1 =CM{iimg};
    Nc=size(CM1,2);
    dist=zeros(winsize,Nc);
    dist1=zeros(winsize,Nc);
    for iwin = 1: winsize
        if iimg<=iwin,
            continue;
        end
        if nargin>3
            dist1(iwin,:) = sqrt(sum(bsxfun(@minus,CM1(1:2,:),newCM(1:2,iimg-iwin)).^2,1));
            dist(iwin,:) = sqrt(sum(bsxfun(@minus,CM1(1:2,:),Cinfo(1:2,iimg-iwin)).^2,1));
        else
            dist(iwin,:) = sqrt(sum(bsxfun(@minus,CM1(1:2,:),Cinfo(1:2,iimg-iwin)).^2,1));
        end
    end
    dist = mean(dist,1);
    dist1 = mean(dist1,1);
    if isempty(dist)
        continue;
    end
    A = zeros(size(MOV,1),size(MOV,2));
    [~, inx]=sort(dist,'ascend');
    if nargin>3        
        [~, inx1]=sort(dist1,'ascend');
        if inx(1) ~= inx1(1),
            continue;
        end
    end
    inx = inx(1);
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
if nargout>2
    MSK = reshape(MSK, size(MOV));
end

