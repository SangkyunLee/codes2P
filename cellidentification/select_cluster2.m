function select_cluster2(cci,Para)
% select_cluster2(cci, Para)
% Usage: (1) mouse draganddrop- remove ROIs in the area from the cci
%        (2) mouse draganddrop with the Shift key - combine ROIS in the area
%        from the cci
%        (3) mouse draganddrop with the Control key - modify/generate ROIs in
%        the area
%        INPUT: meanimg scaled to [0-1], baseimg [mxnx3], mode:'GCamp6' or
%        'OGB' not fully tested, fnsave: savefilename
% Sangkyun Lee - 7/1/2013
% Sangkyun Lee - 9/6/2013; imaging area selection modified.
meanimg = Para.meanimg;
baseimg = Para.baseimg;
mode = Para.mode;
fnsave = Para.fnsave;
figname = Para.figname;
hfig = figure; hist(meanimg(:),100);
uiwait(hfig);
prompt = {'Enter a threshold of imaging area'};
dlg_title = 'Area definition';
num_lines = 1;
def = {'0'};        
answer = inputdlg(prompt,dlg_title,num_lines,def);
        
imagingarea = find(meanimg(:)>str2num(answer{1}));

[dy dx]=size(meanimg);
[dy1 dx1]=size(cci);
[dy2 dx2 ~]=size(baseimg);
if (dy1~=dy2) || (dx1~=dx2)
    error('Dimensions of images should be same');
end

winheight = 700;
winwidth = 1000;
hfig=figure('Position', [400 200 winwidth winheight]);
if exist('figname')
    set(hfig,'Name',figname);
end
left=0.02;  fwidth=0.38; fheight=fwidth*dy/dx;
if fheight>0.5
    fheight=0.5;
    fwidth = fheight*dx1/dy1;
end
fbottom = 1-fheight-0.01;
haxs(1)=axes('Position',[left fbottom fwidth fheight],'Parent',hfig);   
img01=(meanimg-min(meanimg(:)))/(max(meanimg(:))-min(meanimg(:)));
img01 = repmat(img01,[1 1 3]);
hImg=image(img01);
axis equal; xlim([1 dx]); ylim([1 dy]);


left=0.4;  fwidth=0.55; fheight=fwidth*dy1/dx1;
if fheight>0.8
    fheight=0.8;
    fwidth = fheight*dx1/dy1;
end
fbottom = 1-fheight-0.01;
haxs(2)=axes('Position',[left fbottom fwidth fheight],'Parent',hfig);   
hImg2 = overlayImage(baseimg,cci, haxs(2),'Discrete')
axis equal; xlim([1 dx]); ylim([1 dy]);

% display roi number

display_ROInum(haxs(2),cci,cci);

    
    
button_pos0 = [round(0.8*winwidth) round(winheight*(1-((1-fbottom)+0.23))),round(winwidth*0.13) round(winheight*0.06)] 
hbutton0=uicontrol('Style', 'pushbutton', 'String', 'Undo', 'Position', button_pos0); 
set(hbutton0,'Enable','off')    

button_pos1 = [round(0.15*winwidth) round(winheight*(1-((1-fbottom)+0.22))),round(winwidth*0.13) round(winheight*0.06)] 
hbutton1=uicontrol('Style', 'pushbutton', 'String', 'SAVE', 'Position', button_pos1); 
    
button_pos2 = [button_pos1(1)+round(winwidth*0.14) round(winheight*(1-((1-fbottom)+0.22))),round(winwidth*0.13) round(winheight*0.06)] 
hbutton2=uicontrol('Style', 'pushbutton', 'String', 'CLOSE', 'Position', button_pos2); 


handles = guihandles(hfig);
handles.mode = mode;
handles.hfig = hfig;
handles.hbUndo = hbutton0;
handles.hImg2 = hImg2;
handles.winsize=[winwidth winheight];
handles.baseimg = baseimg;
handles.imagingarea = imagingarea;
handles.cci = cci;
handles.newcci = cci;
handles.hax = haxs(2);
handles.fnsave = fnsave;
handles.selectedROIs = [1:length(unique(cci))-1];
handles.removedROIs = {};
handles.Nremoval=0;
handles.bUndo = false;

handles.stackpoint=[];
guidata(hfig,handles);

set(hbutton0,'Callback',{@ButtonDown1});
set(hbutton1, 'Callback', {@savedata}); 
set(hbutton2, 'Callback', {@closef}); 
set(hfig,'WindowButtonDownFcn',{@wbd});

guidata(hfig,handles);
% uiwait(hfig);

function savedata(hObj, event, handles)
    handles = guidata(hObj);
    selectedROIs = handles.selectedROIs;
          
    newcci = handles.newcci;
    choice = questdlg('Would you like to resort ROI order?', ...
        'Re-sort Menu', 'Yes','No','Yes');
    % Handle response
    switch choice
        case 'Yes'
            bresort=1;
        case 'No'
            bresort=0;
    end
    if bresort,    
        listroi = setdiff(unique(newcci(:))',0);
        listloc =[];
        for iroi= 1:length(listroi)
            inxs = find(newcci(:)==listroi(iroi));        
            listloc = [listloc  [listroi(iroi); min(inxs)]];
        end
        [a, order]=sort(listloc(2,:),'ascend');
        for iroi= 1:length(listroi)
            inxs = find(handles.newcci(:)==listloc(1,order(iroi)));
            newcci(inxs)=iroi;
        end
        
    end
    
    
   
    neuropilarea={};
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
        % generate filled cellbodies for the subsequent identification of neuropil
        % area
        filledcells = newcci;
        listroi = setdiff(unique(newcci(:))',0);
        Nroi = length(listroi);

        for iroi = 1: Nroi
            roinum = listroi(iroi);
            inxs = find(filledcells(:)==roinum);    
            [yc,xc] = ind2sub(size(filledcells),inxs);  

            for iy=unique(yc)'
                mm = xc(find(yc==iy));
                filledcells(iy,[min(mm):max(mm)])=roinum;
            end

            for ix=unique(xc)'
                mm = yc(find(xc==ix));
                filledcells([min(mm):max(mm)],ix)=roinum;
            end     
        end
        % identify neuropil area for each cell
        
        prompt = {'Enter pixel resolution(micron/pix)','Enter a radius of neuropil area(micron)'};
        dlg_title = 'Area definition';
        num_lines = 2;
        def = {'1.2','12'};        
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        

        samplingint = str2num(answer{1});
        neuropilradius = round(str2num(answer{2})/samplingint);
        imagingarea=handles.imagingarea;
        nonimagingarea = setdiff([1:prod(size(newcci))],imagingarea);
        for iroi = 1: Nroi
            roinum = listroi(iroi);
            inxs = find(newcci(:)==roinum);
            [yc,xc] = ind2sub(size(newcci),inxs);
            y1 =round((min(yc)+max(yc))/2);
            x1 =round((min(xc)+max(xc))/2);
            y2=y1+[-neuropilradius:neuropilradius];
            y2 = y2(find(y2>0 & y2<size(newcci,1)));
            x2=x1+[-neuropilradius:neuropilradius];    
            x2 = x2(find(x2>0 & x2<size(newcci,2)));

            extractedarea=filledcells(y2,x2);    


            [X Y]=meshgrid([x2-x1],[y2-y1]);
            rho =sqrt(X.^2+Y.^2);
            mask =zeros(size(rho));
            mask(extractedarea(:)==0 & rho(:)<neuropilradius)=1;
            M=zeros(size(newcci));
            M(y2,x2)=mask;
            neuropilarea{roinum}= setdiff(find(M(:)==1),nonimagingarea);
           
        end        
    end
    selectedROIs = setdiff(unique(newcci(:))',0);
    bvalid = true;
    save(handles.fnsave,'selectedROIs','newcci','neuropilarea','bvalid');
    close(handles.hfig);



% guidata(hObject, handles);
function ButtonDown1(hObject, eventdata, handles)
% hObject    handle to axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles = guidata(hObject);
    selectedROIs = handles.selectedROIs;
    removedROIs = handles.removedROIs;
    cci = handles.cci;
%     Obj=get(hObject);
%     hParent=Obj.Parent;


    
    if handles.Nremoval>0  & handles.bUndo
        disp('undo')
        inxROI = removedROIs{handles.Nremoval};
%         removedROIs{handles.Nremoval}=[];
%         selectedROIs = [selectedROIs inxROI];        
        handles.Nremoval = handles.Nremoval-1;
    

        newcci = cci;   
        selectedROIs = setdiff(unique(newcci(:)),0)';
%         removedROIs_mat = cell2mat(removedROIs);
%         for iroi=1:length(removedROIs_mat)
%             inxs = find(newcci(:)== removedROIs_mat(iroi));
%             newcci(inxs)=0;
%         end


        hImg = overlayImage(handles.baseimg,newcci, handles.hax,'Discrete');


        % display roi number
        haxis = gca;
        display_ROInum(haxis,newcci,newcci)

        handles.hImg2 = hImg;
        handles.selectedROIs = selectedROIs;
        handles.removedROIs = removedROIs;
        handles.newcci = newcci;
        handles.bUndo = false;
        set(handles.hbUndo,'Enable','off')
        guidata(hImg, handles);
        
    end
   


function display_ROInum(haxis,cci,newcci,clr)

if nargin<4,
    clr='w';
end    
axes(haxis);
Nroi = length(unique(newcci(:)))-1;
listroi = setdiff(unique(newcci(:)),0);
for iroi= 1:Nroi
    inxs = find(cci(:)==listroi(iroi));
    [yc,xc] = ind2sub(size(newcci),inxs(round(end/2)));
    text(xc,yc,num2str(listroi(iroi)),'Color',clr);
end    

function closef(hObj, event, handles)
    handles = guidata(hObj);
    bvalid = false;
    save(handles.fnsave,'bvalid');
    close(handles.hfig);

    
    

% ---------------------------
function wbd(h,evd)
handles =guidata(h);
Obj = get(h);
handles.fseltype = get(h,'SelectionType');
Xrange=get(handles.hImg2,'Xdata');
Yrange=get(handles.hImg2,'Ydata');
% disp('down')


pos = round(get(Obj.CurrentAxes,'CurrentPoint'));
pos = pos(1,1:2);

if ~(pos(1)<Xrange(1) || pos(1)>Xrange(2) || pos(2)<Yrange(1) || pos(2)>Yrange(2))
%     text(pos(1,1),pos(1,2),'p','Color','y')
    handles.startpoint = pos';


    % set the new values for the WindowButtonMotionFcn and
    % WindowButtonUpFcn
    set(h,'WindowButtonMotionFcn',{@wbm})
    set(h,'WindowButtonUpFcn',{@wbu})
    guidata(h,handles);
end

% ---------------------------
function wbm(h,evd)
% executes while the mouse moves

% disp('motion')

handles =guidata(h);
Xrange=get(handles.hImg2,'Xdata');
Yrange=get(handles.hImg2,'Ydata');



pos = get(get(h,'CurrentAxes'),'CurrentPoint');
pos = round(pos(1,1:2));

if ~(pos(1)<Xrange(1) || pos(1)>Xrange(2) || pos(2)<Yrange(1) || pos(2)>Yrange(2))
% 
    Img = get(handles.hImg2);
    Img =Img.CData;
    handles.hImg2 = image(Img);
    [dy, dx, ~] =size(Img);
    axis equal; xlim([1 dx]); ylim([1 dy]);
    haxis = gca;
    display_ROInum(haxis,handles.newcci,handles.newcci);
    
    
    
    startpoint = handles.startpoint;
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
    
    handles.stackpoint = [pos' handles.stackpoint];
    guidata(h,handles);
end


% ---------------------------
function wbu(h,evd)
% executes when the mouse button is released

disp('up')

handles = guidata(h);
set(h,'WindowButtonMotionFcn','')
set(h,'WindowButtonUpFcn','') 
% guidata(h,handles); 

Xrange=get(handles.hImg2,'Xdata');
Yrange=get(handles.hImg2,'Ydata');

pos = get(get(h,'CurrentAxes'),'CurrentPoint');
pos = round(pos(1,1:2));

if ~(pos(1)<Xrange(1) || pos(1)>Xrange(2) || pos(2)<Yrange(1) || pos(2)>Yrange(2))
    selectedROIs = handles.selectedROIs;
    removedROIs = handles.removedROIs;
    newcci = handles.newcci;
    cci = handles.newcci;
    
    
    startpoint = handles.startpoint;
    endpoint = pos';
    pos_rec=sort([startpoint endpoint],2);
    [X,Y]=meshgrid([pos_rec(1,1):pos_rec(1,2)],[pos_rec(2,1):pos_rec(2,2)]);
    inxs=sub2ind(size(newcci),Y(:),X(:));
    selROIs = setdiff(unique(newcci(inxs)),0);   
    selROIs = selROIs(:)';
    handles.Nremoval = handles.Nremoval+1;                          
    selectedROIs = setdiff(selectedROIs, selROIs);  
    removedROIs {handles.Nremoval} =selROIs;
    
    if strcmp(handles.fseltype,'extend')   
        cellidenties = [unique(newcci(:))' selROIs];
        newclstID =round(rand*200);
        while ~isempty(find(cellidenties == newclstID ))
            newclstID =round(rand*200);
        end
        selectedROIs = [selectedROIs newclstID];
        for iroi=1:length(selROIs)
            inxs = find(newcci(:)== selROIs(iroi));
            newcci(inxs) = newclstID;
        end
    elseif strcmp(handles.fseltype,'alt')
        
        if isfield(handles,'mode') %& strcmp(handles.mode,'GCamp6')
            selinxs=[];
            for iroi=1:length(selROIs)
                subselinxs = find(newcci(:)== selROIs(iroi));
                selinxs = [selinxs; subselinxs(:)];  
                newcci(subselinxs)=0;
            end
            if isempty(selinxs)
                baseimg = handles.baseimg;
                Y1 = startpoint(2)% - 10;
                Y2 = endpoint(2)% + 10;
                X1 = startpoint(1)%-10;
                X2 = endpoint(1)% + 10;
%                 if Y1<1, Y1=1; end
%                 if X1<1, X1=1; end
%                 if X2>size(baseimg,2), X2 = size(baseimg,2); end
%                 if Y2>size(baseimg,1), Y2 = size(baseimg,1); end
                selinxs=sub2ind(size(baseimg),[Y1 Y1 Y2 Y2],[X1 X2 X1 X2]);    
                
%                 figure; image(baseimg([Y1:Y2],[X1:X2],:))
            end
     

            clear global clststruct
            global clststruct
            hnewFig=figure('Position',[100 678 560 420]);
            mansep(hnewFig,handles.baseimg(:,:,1),selinxs,selinxs,selROIs(iroi));
            uiwait(hnewFig);
            if isfield(clststruct,'sel_boundvoxinFOV')            
                Nroi=length(clststruct.sel_boundvoxinFOV);
                for iroi=1:Nroi
                    cellidenties = [unique(newcci(:))' selROIs];
                    inxvoxs = clststruct.sel_boundvoxinFOV{iroi};
                    if ~isempty(inxvoxs)
                        newclstID =round(rand*200);
                        while ~isempty(find(cellidenties == newclstID ))
                            newclstID =round(rand*200);
                        end                
                        newcci(inxvoxs)=newclstID;
                        selectedROIs = [selectedROIs newclstID];
                    end
                end
            end
        end
                
        
        
        
        
        
    else  
        for iroi=1:length(selROIs)
            inxs = find(newcci(:)== selROIs(iroi));
            newcci(inxs)=0;
        end
    end
    
    hImg = overlayImage(handles.baseimg,newcci, handles.hax,'Discrete');


    % display roi number
    haxis = gca;
    display_ROInum(haxis,newcci,newcci)



    handles.hImg2 = hImg;
    handles.newcci = newcci;
    handles.cci = cci;
    handles.selectedROIs = selectedROIs;
    handles.removedROIs = removedROIs;
    handles.bUndo = true;
    set(handles.hbUndo,'Enable','on')
    set(h,'WindowButtonMotionFcn','')
    set(h,'WindowButtonUpFcn','') 
    guidata(h,handles);    
    
end
  

    
