clear all
eye_dir = 'X:\Awake\EYE\eye0129B2_6'
% eye_dir='/media/sdb_WD4T/data_2photon/AWAKE/0822BRLF/11082015_0822BRLF/eye11082015'
recordnum=1
videofn = fullfile(eye_dir,sprintf('%04d.mj2',recordnum));
dlogfn = fullfile(eye_dir,sprintf('%04d.dlog',recordnum));
vlogfn = fullfile(eye_dir,sprintf('%04d.vlog',recordnum));
load(fullfile(eye_dir,sprintf('%04d.mat',recordnum)))
Vreader = VideoReader(videofn)
Vreader.NumberOfFrames



movie = read(Vreader,[1]);
movie_1fr = single(movie(:,:,1));
THR = mean(movie_1fr(:)) +  3*std(movie_1fr(:));

% figure;
% for ii=1:Vreader.NumberOfFrames
%     movie1 = read(Vreader,ii);
%     movie1(movie1(:)>THR) = THR;
%     imagesc(movie1);
%     colormap('gray')
%     pause(0.05);
% end

vfid =fopen(vlogfn,'rt');
vlog=single(fscanf(vfid,'%d, %f, %f\n'));
fclose(vfid)
vlog = reshape(vlog,[3 length(vlog)/3])';


dfid = fopen(dlogfn,'rb');
dlog = fread(dfid,'double');
fclose(dfid);
DAQTriggerTime=dlog(2);
DAQTriggerTime_toc=dlog(3);
dlog =dlog(4:end);
Nch=2
daqdata=reshape(dlog,[Nch length(dlog)/Nch])';

datestr(DAQTriggerTime,'mmmm dd, yyyy HH:MM:SS.FFF AM')
datestr(datenum(data_params.video.InitialTriggerTime),'mmmm dd, yyyy HH:MM:SS.FFF AM')

% datestr(TriggerTime, 'HH:MM:SS.FFF')
% datestr(datenum(data_params.video.InitialTriggerTime), 'HH:MM:SS.FFF')
trigger_offset = str2num(datestr(DAQTriggerTime - datenum(data_params.video.InitialTriggerTime),'SS.FFF'))*1000;



%% eyetracking video analysis
baseimg = squeeze(read(Vreader,1));
baseimg =single(histeq((baseimg)));

par.parname={'Mean'};
par.op1={'>p'};
par.pmaps=baseimg(:,:,1);
par.fnsave='THRtemp.mat';
validate_thresholds(baseimg,par);

Nframe = Vreader.NumberOfFrames;
% 
% load('tmp.mat')
% figure; 
% for ii=2:Nframe
%     a=baseimg;
%     a(PIX{ii}{CInx(ii)})=255;
%     imagesc(a);
%     title(sprintf('%d',ii));
%     pause(0.005);
% end
%     
% for ii=2:Nframe
%  [Y, X]=ind2sub(size(baseimg),PIX{ii}{CInx(ii)});


%% display of pupil identification
hfig = figure;
[dy, dx,~]=size(baseimg);
left=0.1;  

if dy>=dx,
    fheight = 0.7;
    fwidth = fheight*dx/dy;
else
    fwidth = 0.7;
    fheight = fwidth*dy/dx;
end
fbottom = 1 - (fheight+0.02);  
hax=axes('Position',[left fbottom fwidth fheight],'Parent',hfig);  

thr=30;
Cinfo = zeros(3,Nframe);
Frlist = 1 : 10 :Nframe;


j = 1;
for k = 1 : Nframe
    tic
    sMOV1 = read(Vreader,k);
    if mod(k,10)==1,
        if k==1,
            sMOV = uint8(zeros(size(sMOV1,1),size(sMOV1,2),length(Frlist)));
        end
        sMOV(:,:,j) = histeq(sMOV1);
        j = j+1;
    end
    ti = toc;
    fprintf('T(%d):%f\n',k,ti);
end


for iimg = 1 : size(baseimg1,3)

    Imgi = baseimg1(:,:,iimg);
    fmap = zeros(size(Imgi));
    fmap(Imgi<thr)=1;
    if iimg==1
        fmap1=fmap;
        [L1, num] = bwlabel(double(fmap1),8);
        hfig1 = figure; imagesc(L1);
        uiwait(hfig1);
        prompt = {'Enter pupil cluster ID'};
        dlg_title = 'PUPIL ID';
        num_lines = 1;
        def = {'1'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        ID = str2double(answer{1});
        L1(L1(:)~=ID)=0;
        Cpix = find(L1(:)==ID);
        [X0 Y0]=cmass(L1);
        CM0=[X0 Y0 length(Cpix)];
        
        
    else
        [Li num] = bwlabel(double(fmap),8);
        for ic=1:num
            cpix = find(Li(:)==ic);
            if length(cpix)<100
                Li(cpix)=0;
            end
        end
        CID = setdiff(unique(Li(:)),0);
        Nc = length(CID);
        CM=zeros(3,Nc);
        for ic = 1 : Nc
            Lx = zeros(size(Li));
            Cpix1 = find(Li(:)==CID(ic));
            Lx(Cpix1) = ic;
            [X0 Y0]=cmass(Lx);
            CM(:,ic)=[X0 Y0 length(Cpix1)];            
        end
        winsize=20;
        dist=zeros(winsize,Nc);
        for iwin = 1: winsize
            if iimg<=iwin,
                continue;
            end
            dist(iwin,:,iimg) = sqrt(sum(bsxfun(@minus,CM(1:2,:),Cinfo(1:2,iimg-iwin)).^2,1));
        end
        dist = mean(dist,1);
        %ds =abs(bsxfun(@minus,CM(3,:),Cinfo(3,iimg-1)));
        [~, inx]=min(dist);
        CM0=CM(:,inx);
        
        L1 = zeros(size(Li));
        L1(Li(:)==CID(inx)) =1; 
        
    end
    Cinfo(:,iimg) = CM0;

    baseimg1 = repmat(baseimg1, [1 1 3]);
    overlayImage(baseimg1,L1, hax,'Continuous');
    title(sprintf('%d',iimg));
    pause(0.002);    
end

Xr=floor([min(Cinfo(1,1:iimg-1))-100 max(Cinfo(1,1:iimg-1))+100]);
Yr=floor([min(Cinfo(2,1:iimg-1))-100 max(Cinfo(2,1:iimg-1))+50]);
figure; 
H=fspecial('gaussian',5,5);
for iimg = 1 : Nframe
    baseimg1 = read(Vreader, iimg);
    baseimg1 = single(histeq(baseimg1));
    tmpImg= baseimg1(Yr(1):Yr(2),Xr(1):Xr(2));
    if iimg ==1,
        PImg=zeros(size(tmpImg,1),size(tmpImg,2),Nframe);
        fImg=zeros(size(tmpImg,1),size(tmpImg,2),Nframe);
    end
    PImg(:,:,iimg)=tmpImg;
    fImg(:,:,iimg)=imfilter(tmpImg,H,'same');
    subplot(121);
    imagesc(histeq(PImg(:,:,iimg)));axis image;
    subplot(122);
    imagesc(histeq(fImg(:,:,iimg)));axis image;
    pause(0.005);
    
end

figure;
for iimg=1:Nframe
    subplot(121);
    imagesc(histeq(PImg(:,:,iimg)));axis image;
    subplot(122);
    imagesc(histeq(fImg(:,:,iimg)));axis image;
    pause(0.005);
end

    
%% 
% Thr=12
% stframe = 1001
% Nframe = 1000;
% fmapall=single(zeros(size(baseimg,1),size(baseimg,2),Nframe));
% movie = read(Vreader,[stframe stframe+Nframe-1]);
% movie = squeeze(movie(:,:,1,:));
% fmap = single(zeros(size(movie)));
% fmap(movie(:)<Thr)=1;
% 
% [L,NUM] = spm_bwlabel(double(fmap),26);
% for ii=1:NUM
%     aa=find(L(:)==ii);
%     if length(aa)<1000,
%         L(aa)=0;
%     end
% end
% L(L(:)~=63)=0;
%     
% for ii=1:Nframe
% imagesc(L(:,:,ii)); colorbar; caxis([0 3])
% title(num2str(ii),'FontSize',20);
% pause(0.01);
% end
% 
% L(L(:)>3)=3;
% figure; hax=subplot(1,1,1);
% for iimg = 1 : Nframe
%     baseimg1 = movie(:,:,iimg);
%     baseimg1(baseimg1(:)>50)=50;
%     baseimg1 = double(baseimg1(:,:,1));
% 
%     baseimg1 = repmat(baseimg1, [1 1 3]);
%     overlayImage(baseimg1,L(:,:,iimg), hax,'Discrete',[1 3]);
%     title(num2str(iimg),'FontSize',20);
%     caxis([0 3]);
%     pause(0.01);    
% end
% 
