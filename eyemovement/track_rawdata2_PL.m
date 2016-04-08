function [out, MSK]= track_rawdata2_PL(MOV,opts, newCM,htwin, frames)
% function track_rawdata(MOV,opts)
% opts.winsize=10;
% opts.thr =30
% opts.bdisp =false;


winsize = opts.winsize; 
thr = opts.thr;
if nargin<5
    htwin=0;
end



if nargout>2
    MSK = uint8(zeros(size(MOV,1)*size(MOV,2),size(MOV,3)));    
end
%% create clusters
Nframe = size(MOV,3);

CM = cell(Nframe,1);
PIX = cell(Nframe,1);



Nsel = length(frames);
CM_sel = cell(Nsel,1);
PIX_sel = cell(Nsel,1);
MSK_sel = uint8(zeros(size(MOV,1)*size(MOV,2),Nsel));  

parfor iimg0 = 1: Nsel
    MOV1 =MOV;
 
    iimg = frames(iimg0);    


    MSK1 = zeros(size(MOV1,1),size(MOV1,2));
    fmap = zeros(size(MOV1,1),size(MOV1,2));
    if (iimg-htwin)<1 || (iimg+htwin)>Nframe
        a=MOV1(:,:,iimg);
        fmap(a(:)<thr)=1;
    else
        a=mean(MOV1(:,:,iimg-htwin:iimg+htwin),3);
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

    PIX_sel{iimg0}=PIX0;
    CM_sel{iimg0}=CM0;
    MSK_sel(:,iimg0)=uint8(MSK1(:));
    
end
PIX(frames) = PIX_sel(:);
CM(frames) = CM_sel(:);
MSK(:,frames) = MSK_sel;
clear PIX_sel CM_sel MSK_sel;
pause(1);

%% find the closest cluster as the pupil  
sPIX = cell(Nframe,1);
% CInx = zeros(Nframe,1);
Cinfo = zeros(2,Nframe);
dists = zeros(Nframe,1);

if opts.bdisp
    figure; 
end
for iimg = frames  

    CM1 =CM{iimg};
    Nc=size(CM1,2);
    dist=zeros(winsize,Nc);    
    for iwin = 1: winsize
        if iimg <= iwin,
            continue;
        end
        dist(iwin,:) = sqrt(sum(bsxfun(@minus,CM1(1:2,:),newCM(1:2,iimg-iwin)).^2,1));            
    end

    dist = mean(dist,1);
    
    A = zeros(size(MOV,1),size(MOV,2));
    [~, inx]=min(dist);
    
    if ~isempty(inx)
        dists(iimg)=dist(inx);
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
        
    end

    
        
    
    if opts.bdisp
        imagesc(A);
        title(num2str(iimg));
        pause(0.002);
    end
end

out.sPIX = sPIX;
out.CM = Cinfo;
out.validframes = frames;
if nargout>2
    MSK = reshape(MSK, size(MOV));
end

