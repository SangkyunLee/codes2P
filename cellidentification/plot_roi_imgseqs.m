function plot_roi_imgseqs(mainpath, fnpattern, numdata,ROI)
% mainpath = 'D:\data_2photon\140103 BR\RF_140103BR_S1_20140206\data'
% fnpattern = 'TSeries-02062014-1233-001';
% numdata = [1:22];

fn_MC=cell(length(numdata),10);
data_subpath = cell(1,max(numdata));
data_Nseg=zeros(1, 10);
for idata =1 : length(numdata)  
    data_subpath{idata}= sprintf([fnpattern '-%03d'],numdata(idata));
    flist = dir(fullfile(mainpath,data_subpath{idata}));
    iseg = 1;
    for ifl= 1 : length(flist)
        if ~flist(ifl).isdir
            if strcmp(flist(ifl).name(end-3:end),'.tif') && isempty(strfind(flist(ifl).name,'Voltage'))            
                fn1= sprintf(['MC_' data_subpath{idata} '_Cycle00001_Ch2-%02d.zip'],iseg);                
                fn2 = fullfile(mainpath,data_subpath{idata},fn1);
                if exist(fn2,'file')
                    fn_MC{idata,iseg} =fn1;
                    iseg = iseg +1;
                end
            end
        end
    end
    data_Nseg(idata)=iseg-1;
            
end
%------------------------------

% M.data.ROI = newROI;
for idata = 1 : length(numdata)  
    for iseg = 1: data_Nseg(idata)
        fnfig=sprintf('initialROI%d-%02d.jpg',numdata(idata),iseg);
        subpath=sprintf([fnpattern '-%03d'],numdata(idata));
        fn=sprintf(['AVG_MC_' fnpattern '-%03d_Cycle00001_Ch2-%02d.tif'],numdata(idata),iseg);
        finfo=fullfile(mainpath,subpath,fn);
        if ~exist(finfo,'file')
            fn=sprintf(['MC_' fnpattern '-%03d_Cycle00001_Ch2-%02d.tif'],numdata(idata),iseg);
            finfo=fullfile(mainpath,subpath,fn);
            bAVG = false;
        else
            bAVG = true;
        end
            
        
        Img =loadTseries(finfo);
        Img = single(Img);
        img1 = mean(Img,3);
        if ~bAVG
            avgfn = sprintf(['AVG_MC_' fnpattern '-%03d_Cycle00001_Ch2-%02d.tif'],numdata(idata),iseg);
            fullfn_avgMC = fullfile(mainpath,subpath,avgfn);
            imwrite(uint16(img1),fullfn_avgMC,'Compression','none');            
        end

        %[ysize xsize]=size(img1);
        mimg1 = mean(img1(:)); stdimg1 = std(img1(:));
        thr = round(mimg1+5*stdimg1);
        img1(img1(:)>thr) = thr;
        baseimg = img1;


        [m, n, ~] = size(baseimg);
        mask =zeros(m*n,1);
        nROI = length(ROI);
        for ir=1:nROI            
            inx = ROI(ir).boundary_list;
            mask(inx)=1;
        end

        %mask = reshape(mask, [m n]);


        hfig=figure;                 
        imagesc(baseimg); colormap('gray');
        set(hfig,'Position',[0 0 1024 1024])
        axis image
        if exist('./ROI','dir')==0,
            mkdir('ROI')
        end
        set(gca,'position',[0 0 1 1],'units','normalized')
        print('-dtiff',['./ROI/0-' fnfig]);
        pause(1);
%         close(hfig)
%         opts.alphalev=1
%         [~, hfig] = overlayImage2(M.data.img1,mask,opts);
%         
        nROI = length(ROI);        
        for iroi= 1:nROI
            if ~isempty(ROI(iroi).centerPos)
                yc = ROI(iroi).centerPos(1);    
                xc = ROI(iroi).centerPos(2);        
                text(xc,yc,ROI(iroi).name,'Color','r','FontSize',5);    
            end
        end
        axis image; 
        set(hfig,'Position',[0 0 1024 1024])
        if exist('./ROI','dir')==0,
            mkdir('ROI')
        end


        
        set(gca,'position',[0 0 1 1],'units','normalized')
        print('-dtiff',['./ROI/' fnfig]);
        %saveas(hfig,['.\\ROI\\' fnfig],'tif')
        
        pause(1);
        close(hfig)
    end
end