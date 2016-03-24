function M=ROIoverlap(img1,ROI_list)
% function M=ROIoverlap(img1,ROI_list)
% INPUT: 
%   img1: meanimage
%   ROI_list: structure containing roi info
%              
% 2013-10-23, written by Sangkyun Lee
% 2013-12-04, modified by Sangkyun Lee
% 2014-01-22, Add delete function, modified by Sangkyun Lee
% 2015-09-28, bug fixed by Sangkyun Lee
[ysize xsize]=size(img1);
M.data.img1 = img1;
ROI_list = formulateROI(ROI_list,[ysize xsize]);
M.data.ROI = ROI_list;
M.data.newROI = ROI_list;
M.data.btext=true;

xWinsize = 700;
yWinsize = 500;


%% set Diplay Range
hfig = figure('Name','Histogram of image pixels'); hist(M.data.img1(:),100);
prompt = {'Enter display range'};
dlg_title = 'Display Range';
num_lines = 1;
maxpix = max(M.data.img1(:));
def_str = sprintf('[100 %s]',num2str(maxpix));
def = {def_str};        
answer = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answer)
    answer=def;
end
eval(['M.ui.DispRange=' answer{1}]);
if exist('hfig','var')
    close(hfig);
end

%% create figure;
figName = sprintf('Overlap ROIs ');
M.fig = figure('Color', 'w', 'Name', figName, ...
			   'Units', 'pixels', 'Position', [100 100 xWinsize yWinsize]);

set(M.fig,'WindowButtonDownFcn',{@wbd});
%% create axes and uis	   
M.ui.imgAxes = subplot('Position', [.05 .05 .7 .9]);
M.ui.hgp1 = uibuttongroup('Units','normalized','Position',[.78 .68 .20 .2],'parent',M.fig,'SelectionChangeFcn',@selchangeCallback);
M.ui.rad1 = uicontrol('Style','Radio','Units','normalized','FontSize',11,...
    'String','Out-bound','position',[.1 .7 .8 .3],...
    'parent',M.ui.hgp1,'HandleVisibility','off');
M.ui.rad2 = uicontrol('Style','Radio','Units','normalized','FontSize',11,...
    'String','In-bound','position',[.1 .4 .8 .3],...
    'parent',M.ui.hgp1,'HandleVisibility','off');

M.ui.rad3 = uicontrol('Style','Radio','Units','normalized','FontSize',11,...
    'String','None','position',[.1 .1 .8 .3],...
    'parent',M.ui.hgp1,'HandleVisibility','off');


M.ui.xpos = mrvSlider([.78 .58 .2 .08], 'Move_X', ...
    'Range', [-xsize xsize], 'IntFlag', 1, 'Value', 0, ...
    'FontSize', 11, 'Color', 'w','Callback','updateROIdisp');

M.ui.ypos = mrvSlider([.78 .43 .2 .08], 'Move_Y', ...
    'Range', [-ysize ysize], 'IntFlag', 1, 'Value', 0, ...
    'FontSize', 11, 'Color', 'w','Callback','updateROIdisp');

M.ui.setposnest = uicontrol('Style','pushbutton','Units','normalized',...
    'position',[.78 .35 .15 .07],...
    'String','Set Pos & Re-estimate','Callback',@setposnestCallback);

M.ui.export = uicontrol('Style','pushbutton','Units','normalized',...
    'position',[.78 .25 .15 .07],...
    'String','EXPORT','Callback',@exportCallback,'Enable','on');

M = gen_ROIimg(M);


updateROIdisp(M);
end

%%
function exportCallback(varargin)
    M = get(gcf, 'UserData');
    ROI = M.data.newROI;
    choice = questdlg('Would you like to create neuropil area?', ...
        'Neuropil area', 'Yes','No','Yes');
    % Handle response
    switch choice
        case 'Yes'
            bneuropil=1;
        case 'No'
            bneuropil=0;
    end
    
    if bneuropil,
        
        
        hfig = figure('Name','Histogram of image pixels'); hist(M.data.img1(:),100);
        prompt = {'Enter pixel resolution(micron/pix)','Enter a radius of neuropil area(micron)','Set threshold for background image'};
        dlg_title = 'Area definition';
        num_lines = 3;
        hthr = round(mean(M.data.img1(:))+std(M.data.img1(:))*6);
        %def = {'0.518','8, 20',num2str(hthr)};        
        def = {'0.9','7, 20',num2str(hthr)};        
        answer = inputdlg(prompt,dlg_title,num_lines,def);        
        close(hfig);
        
        
        samplingint = str2double(answer{1});
        pixelres = samplingint;
        
        remain = answer{2};
        delim = strfind(remain,',');
        if ~isempty(delim)
            delim = strfind(remain,' ');
        end            
        
        if ~isempty(delim)
            boundary=[];
            ii = 1;
            while true
                [str, remain] =  strtok(remain,',')
                if isempty(str), break; end
                boundary(ii) = str2double(str);
                ii = ii + 1;
            end        
            neuropilradiusinpix = round(boundary(2)/samplingint);
            outerbound = round(boundary(2)/samplingint);
            innerbound = round(boundary(1)/samplingint);
        else
            neuropilradiusinpix = round(str2double(answer{2})/samplingint);
            outerbound = neuropilradiusinpix;
            innerbound = 0;
        end
        
        nROI = length(ROI);
        [m, n] = size(M.data.img1);
        CellareaImg = zeros(m,n);
        listroi = zeros(1,nROI);        
        outercellinpix = innerbound;
        
        
        % to remove any potential contamination from not the most adajecent
        % cellbodies, define cellboundary with the innerbound of the
        % neuropil patch
        for iroi = 1: nROI
            CellareaImg(ROI(iroi).fill_list)=iroi;
            listroi(iroi) = iroi;
            if ~isempty(ROI(iroi).centerPos)
                y1 = round(ROI(iroi).centerPos(1));
                x1 = round(ROI(iroi).centerPos(2));
                y2=y1+(-outercellinpix:outercellinpix);
                y2 = y2(y2>0 & y2<m);
                x2=x1+(-outercellinpix:outercellinpix);    
                x2 = x2(x2>0 & x2<n);

                extractedarea = CellareaImg(y2,x2);    
                [X Y] = meshgrid(x2-x1,y2-y1);
                rho = sqrt(X.^2+Y.^2);
                mask = zeros(size(rho));
                mask(extractedarea(:)==0 & rho(:) < outercellinpix)=1;
                tmpImg = zeros(m,n);
                tmpImg(y2,x2) = mask; 

                CellareaImg = CellareaImg+tmpImg;
            end
            
        end
        
        %--------- define a stable imaging area
        par.parname={'Mean'};
        par.op1={'<p'};
        IMG= M.data.img1;
        hthr = str2double(answer{3});
        IMG(IMG(:)>hthr)=hthr;
        par.pmaps=M.data.img1;
        par.fnsave='THRtemp.mat';
        validate_thresholds(IMG,par);
        load(par.fnsave);        
        imagingarea = THR.fmap;
        ROI(1).thrs = THR.thrs;
        for iroi = 1: nROI
            if ~isempty(ROI(iroi).centerPos)
                y1 = round(ROI(iroi).centerPos(1));
                x1 = round(ROI(iroi).centerPos(2));

                y2=y1+(-neuropilradiusinpix:neuropilradiusinpix);
                y2 = y2(y2>0 & y2<m);
                x2=x1+(-neuropilradiusinpix:neuropilradiusinpix);    
                x2 = x2(x2>0 & x2<n);


                extractedarea=CellareaImg(y2,x2);    


                [X Y]=meshgrid(x2-x1,y2-y1);
                rho =sqrt(X.^2+Y.^2);
                mask =zeros(size(rho));
                mask(extractedarea(:)==0 & rho(:) <= outerbound & rho(:)> innerbound)=1;
                tmpImg2=zeros(m,n);
                tmpImg2(y2,x2)=mask;
                ROI(iroi).neuropilarea = find(tmpImg2(:)& imagingarea(:));
            else
                ROI(iroi).neuropilarea = [];
            end
           
        end 
    else
        pixelres = [];
        neuropilradiusinpix = [];
    end    
%     choice = questdlg('Would you like to sort ROIs?', ...
%         'Sorting ROIs', 'Yes','No','No');
%     
%     switch choice
%         case 'Yes'
%             bsort=1;
%         case 'No'
%             bsort=0;
%     end
%     if bsort
%         [m, n] = size(M.data.img1);
%         nROI = length(ROI);
%         Pos = zeros(2, nROI);
%         for iroi = 1:nROI
%             Pos(:,iroi) = ROI(iroi).centerPos(:);
%         end
%         inx = sub2ind([m,n],round(Pos(1,:)),round(Pos(2,:)));
%         [mv newinx]=sort(inx,'ascend');
% 
%         clear newROI;
%         for iroi=1:nROI    
%             newROI(iroi) = ROI(newinx(iroi));
%             newROI(iroi).name = num2str(iroi);
%         end
%         clear ROI;
%         ROI = newROI;
%     
%     end

    
    assignin('base','newROIs',ROI);
    assignin('base','pixelres',pixelres);
    assignin('base','neuropilradiusinpix',[innerbound outerbound]);
    set(M.ui.export,'Enable','off');
    close(M.fig);
end


%% set the new position and re-estimate rois
function setposnestCallback(varargin)
    M = get(gcf, 'UserData');

    set(M.ui.setposnest,'Enable','off');
    set(M.ui.export,'Enable','off');
    set(M.fig,'UserData',M);
    drawnow;
    xmove = get(M.ui.xpos.sliderHandle, 'Value');
    ymove = get(M.ui.ypos.sliderHandle, 'Value');
    ROI = M.data.newROI;
    [m, n]=size(M.data.img1);
    nROI = length(ROI);
%     newROI=[];
    for iroi = 1:nROI
        iroi
        
        oldPos = ROI(iroi).centerPos;
        if ~isempty(oldPos)
            currPos = oldPos + [ymove xmove];        

            [Y X]=ind2sub([m n],ROI(iroi).fill_list);
            dist =sqrt((Y-oldPos(1)).^2 + (X-oldPos(2)).^2);
            cell_radius = max(dist);
            pixel_dimension=[1 1];
            clear tmpROI;
            if isfield(ROI(iroi),'type') 
                if strcmp(ROI(iroi).type,'Ring')
                    pixel_list=donut_roi2(M.data.img1,currPos,cell_radius,pixel_dimension,0);
                    tmpROI.type = 'Ring';
                elseif strcmp(ROI(iroi).type,'Disk')
                    pixel_list=disk_roi2(M.data.img1,currPos,cell_radius,pixel_dimension,0);
                    tmpROI.type = 'Disk';
                else
                    error('unspecified type');
                end
            else
                error('ROI type is required');
            end

            if ~isempty(pixel_list)
                tmpROI.pixel_list=pixel_list;            
            else
                tmpROI.pixel_list = ROI(iroi).pixel_list;            
                warning([mfilename ':setposnetCallback'],'estimation failed in ROI %d',(iroi));

            end
            
        
            
        end
        tmpROI = formulateROI(tmpROI, [m n]);            
        tmpROI.name = ROI(iroi).name;    
        tmpROI.type = ROI(iroi).type;
        
        if iroi==1, 
            newROI = tmpROI;
        else
            newROI(iroi) =  tmpROI;
        end
    end
    M.data.newROI = newROI;    
    M=gen_ROIimg(M);
    set(M.ui.xpos.sliderHandle, 'Value',0);
    set(M.ui.ypos.sliderHandle, 'Value',0);
    set(M.ui.xpos.editHandle, 'String',0);
    set(M.ui.ypos.editHandle, 'String',0);
    updateROIdisp(M);
    drawnow();
    set(M.ui.export,'Enable','on');
    set(M.ui.setposnest,'Enable','on');
    set(M.fig,'UserData',M);
    
end


%% 
function selchangeCallback(varargin)
M = get(gcf, 'UserData'); 
M=gen_ROIimg(M);
    updateROIdisp(M);
end
%% 
function M=gen_ROIimg(M)
%     if ~exist('M','var') || isempty(M), M = get(gcf, 'UserData'); end
    [m, n, ~] = size(M.data.img1);
    mask =zeros(m*n,1);
    nROI = length(M.data.newROI);
    switch get(get(M.ui.hgp1,'SelectedObject'),'String')
        case 'Out-bound'        
            for ir=1:nROI            
                inx = M.data.newROI(ir).boundary_list;
                mask(inx)=1;
            end
            M.data.btext=true;
        case 'In-bound'         
            for ir=1:nROI
                inx = M.data.newROI(ir).pixel_list;
                mask(inx)=1;
            end
            M.data.btext=true;
        case 'None'
            M.data.btext=false;
            %nothing
        otherwise
            error('unexpected type');
    end
    mask = reshape(mask, [m n]);
    M.data.roiimg = mask(:,:);     
    set(M.fig,'UserData',M);
end
%%
function out = formulateROI(ROI,Msize)
    nROI = length(ROI);
   
    for i=1:nROI
        tempROI.pixel_list=ROI(i).pixel_list;
        if isempty(tempROI.pixel_list)
           tempROI.boundary_list=[];
           tempROI.centerPos=[];
           if isfield(ROI(i),'name')&& (~isempty(ROI(i).name))
               tempROI.name=ROI(i).name;
           else
               tempROI.name=[];
           end

            if isfield(ROI(i),'type')
                tempROI.type=[];        
            end
            if isfield(ROI(i),'fmean')
                tempROI.fmean=[];        
            end            
            tempROI.fill_list=[];
            if isfield(ROI(i),'fmean_fill')
                tempROI.fmean_fill=[];        
            end
           
        else
            ROImap=false(Msize(1),Msize(2));
            ROImap(ROI(i).pixel_list)=1;

    %         if isfield(ROI(i),'boundary_list')&& (~isempty(ROI(i).boundary_list))
    %             tempROI.boundary_list=ROI(i).boundary_list;
    %         else
                tempROI.boundary_list=find(bwperim(ROImap,4));
    %         end


            [I1,J1]=ind2sub([Msize(1) Msize(2)],tempROI.pixel_list);
            tempROI.centerPos=[mean(I1),mean(J1)];


    %         if isfield(ROI(i),'seedPos')&& (~isempty(ROI(i).seedPos))
    %             tempROI.seedPos=ROI(i).seedPos;
    %         else
    %             tempROI.seedPos=tempROI.centerPos;
    %         end

            if isfield(ROI(i),'name')&& (~isempty(ROI(i).name))
                tempROI.name=ROI(i).name;
            else
                tempROI.name='ROI';
            end

             if isfield(ROI(i),'type')
                tempROI.type=ROI(i).type;        
             end

             if isfield(ROI(i),'fmean')
                tempROI.fmean=ROI(i).fmean;        
             end


            ROImap_fill=imfill(ROImap,'holes');
            tempROI.fill_list=find(ROImap_fill);

            if isfield(ROI(i),'fmean_fill')
                tempROI.fmean_fill=ROI(i).fmean_fill;        
            end
        end
        
        out(i) = tempROI;

    end
end







% ---------------------------
function wbd(h,~)
    M = get(h,'UserData');    
    [m, n] = size(M.data.img1);
    % disp('down')


    pos = round(get(M.ui.imgAxes,'CurrentPoint'));
    pos = pos(1,1:2);
    
    if ~(pos(1)<1 || pos(1)>n || pos(2)<1 || pos(2)>m)
        M.ui.mouse.startpoint = pos';
        M.ui.mouse.stackpoint = [];
        % set the new values for the WindowButtonMotionFcn and
        % WindowButtonUpFcn
        set(h,'WindowButtonMotionFcn',{@wbm})
        set(h,'WindowButtonUpFcn',{@wbu})
        set(h,'UserData',M);
    end
end

% ---------------------------
function wbm(h,~)
    % executes while the mouse moves

    M = get(h,'UserData'); 
    [m, n] = size(M.data.img1);




    pos = round(get(M.ui.imgAxes,'CurrentPoint'));
    pos = round(pos(1,1:2));
    if ~(pos(1)<1 || pos(1)>n || pos(2)<1 || pos(2)>m)
        opts.hAx = M.ui.imgAxes;
        hImg = overlayImage2(M.data.img1, M.data.roiimg, opts);

        startpoint = M.ui.mouse.startpoint;
        endpoint = pos';
        Xpos = [startpoint(1) endpoint(1)];
        Ypos = [startpoint(2) endpoint(2)];

        lineX = [Xpos(1) Xpos(1)];
        lineY = [Ypos(1) Ypos(2)];
        line(lineX,lineY,'Color','y');
        lineX = [Xpos(1) Xpos(2)];
        lineY = [Ypos(2) Ypos(2)];
        line(lineX,lineY,'Color','y');
        lineX = [Xpos(2) Xpos(2)];
        lineY = [Ypos(1) Ypos(2)];
        line(lineX,lineY,'Color','y');
        lineX = [Xpos(1) Xpos(2)];
        lineY = [Ypos(1) Ypos(1)];
        line(lineX,lineY,'Color','y');

        M.ui.mouse.stackpoint = [pos' M.ui.mouse.stackpoint];
        set(h,'UserData',M);
    end
end


% ---------------------------
function wbu(h,~)
% executes when the mouse button is released

    disp('up')

    M = get(h,'UserData');   
    set(h,'WindowButtonMotionFcn','')
    set(h,'WindowButtonUpFcn','') 
    
    ROI = M.data.newROI;
    [m, n] = size(M.data.img1);

    pos = round(get(M.ui.imgAxes,'CurrentPoint'));
    pos = round(pos(1,1:2));

    if ~(pos(1)<1 || pos(1)>n || pos(2)<1 || pos(2)>m)


        startpoint = M.ui.mouse.startpoint;
        endpoint = pos';
        pos_rec=sort([startpoint endpoint],2);
        nROI = length(M.data.newROI);
        selROIs =zeros(1,nROI);
        for iroi = 1:nROI
            if ~isempty(M.data.newROI(iroi).centerPos)
                Yc = M.data.newROI(iroi).centerPos(1);
                Xc = M.data.newROI(iroi).centerPos(2);
                if pos_rec(1,1)<=Xc && pos_rec(1,2)>=Xc && pos_rec(2,1)<=Yc && pos_rec(2,2)>=Yc
                    selROIs(iroi)=1;
                end
            end
        end
        selROIs = find(selROIs);        
        strroi = num2str(selROIs);
        if isempty(selROIs)
           
            prompt = {'Enter a new ROI number','Ring(0) or Disk(1)'};
            dlg_title = 'Select a ROI type';
            num_lines = 2;
            def = {num2str(length(ROI)+1),'0'};        
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            if ~isempty(answer)
                roinum = str2double(answer{1});
                if roinum ==0,
                    return;
                end
                roitype = str2double(answer{2});

                cPos0 = mean(pos_rec,2);
                cPos = [cPos0(2); cPos0(1)];
                bnew = true;
                bcancel = false;
            else
                
                bcancel = true;
            end
        else
            prompt = {['Select ONE from ROIs(' strroi ')'],'Ring(0) or Disk(1)/ Delete(2)'};
            dlg_title = 'Select a ROI type';
            num_lines = 2;
            def = {num2str(selROIs(1)),'0'};        
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            if ~isempty(answer)
                roinum = str2double(answer{1});
                roitype = str2double(answer{2});

                cPos0 = mean(pos_rec,2);
                cPos = [cPos0(2); cPos0(1)];
                bnew = false;
                bcancel = false;
            else
                bcancel = true;
            end
        end
        if ~bcancel
            if roitype==2,
                newROI = M.data.newROI;
                for iROI=1:length(ROI)
                    if strcmp(ROI(iROI).name,num2str(roinum))
                        newROI(iROI).pixel_list=[];
                        newROI(iROI).boundary_list=[];
                        newROI(iROI).centerPos=[];
                        newROI(iROI).fill_list=[];
                        if isfield(newROI(iROI),'neuropilarea')
                            newROI(iROI).neuropilarea=[];
                        end
                        
                        
                    end
                end
                
                M.data.newROI = newROI;
            else   
                [X,Y]=meshgrid(pos_rec(1,1):pos_rec(1,2),pos_rec(2,1):pos_rec(2,2));
                dist =sqrt((Y-cPos(1)).^2 + (X-cPos(2)).^2);
                cell_radius = max(dist(:));
                pixel_dimension=[1 1];
                clear tmpROI
                if roitype==1
                    pixel_list=disk_roi2(M.data.img1,cPos,cell_radius,pixel_dimension,0);
                    tmpROI.type = 'Disk';            
                elseif roitype==0
                    pixel_list=donut_roi2(M.data.img1,cPos,cell_radius,pixel_dimension,0);
                    tmpROI.type = 'Ring';            
                else
                    pixel_list =[];
                end

                if ~isempty(pixel_list)
                    tmpROI.pixel_list=pixel_list;            
                else        
                    tmpROI.pixel_list = ROI(roinum).pixel_list;            
                    msg = sprintf('estimation failed in ROI: %d',iroi);
                    warning(msg);
                end
                if ~bnew
                    tmpROI.name = ROI(roinum).name;                
                else
                    tmpROI.name = num2str(roinum);
                end
                if isfield(M.data.newROI(1),'fmean')
                    tmpROI.fmean =[];
                end
                if isfield(M.data.newROI(1),'fmean_fill')
                    tmpROI.fmean_fill =[];
                end
                tmpROI = formulateROI(tmpROI, [m n]);            
                if (length(M.data.newROI)<roinum && isfield(M.data.newROI(roinum-1),'neuropilarea'))...
                    || (length(M.data.newROI)>=roinum && isfield(M.data.newROI(roinum),'neuropilarea'))
                    tmpROI.neuropilarea=[];
                end
                M.data.newROI(roinum)= tmpROI;
            end
        end
        M=gen_ROIimg(M);
        updateROIdisp(M);

        
        set(h,'WindowButtonMotionFcn','')
        set(h,'WindowButtonUpFcn','') 
        set(h,'UserData',M);

    end
end
