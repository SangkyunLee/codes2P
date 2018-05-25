function select_regions(fmap,Para)

meanimg = Para.meanimg;
baseimg = Para.baseimg;
mode = Para.mode;
fnsave = Para.fnsave;
figname = Para.figname;

[dy dx]=size(meanimg);
[dy1 dx1]=size(fmap);
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
left=0.02; bottom=0.45; fwidth=0.38; fheight=fwidth*dy/dx;
haxs(1)=axes('Position',[left bottom fwidth fheight],'Parent',hfig);   
hImg1 = overlayImage(baseimg,fmap, haxs(1),'Continuous')
axis equal; xlim([1 dx]); ylim([1 dy]);

left=0.3;  bottom=0.22; fwidth=0.65; fheight=fwidth*dy1/dx1;
haxs(2)=axes('Position',[left bottom fwidth fheight],'Parent',hfig);   
hImg2 = overlayImage(baseimg,fmap, haxs(2),'Continuous')
axis equal; xlim([1 dx]); ylim([1 dy]);

    
    
button_pos0 = [round(0.8*winwidth) round(winheight*(1-(fheight+0.23))),round(winwidth*0.13) round(winheight*0.06)] 
hbutton0=uicontrol('Style', 'pushbutton', 'String', 'Undo', 'Position', button_pos0); 
set(hbutton0,'Enable','off')    

button_pos1 = [round(0.24*winwidth) round(winheight*(1-(fheight+0.22))),round(winwidth*0.13) round(winheight*0.06)] 
hbutton1=uicontrol('Style', 'pushbutton', 'String', 'SAVE', 'Position', button_pos1); 
    
button_pos2 = [button_pos1(1)+round(winwidth*0.14) round(winheight*(1-(fheight+0.22))),round(winwidth*0.13) round(winheight*0.06)] 
hbutton2=uicontrol('Style', 'pushbutton', 'String', 'CLOSE', 'Position', button_pos2); 


handles = guihandles(hfig);
handles.mode = mode;
handles.hfig = hfig;
handles.hbUndo = hbutton0;
handles.hImg2 = hImg2;
handles.winsize=[winwidth winheight];
handles.baseimg = baseimg;
handles.fmap = fmap;
handles.cci = zeros(size(fmap));
handles.hax = haxs(2);
handles.fnsave = fnsave;
handles.selectedROIs = [];
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
          
    newcci = handles.cci;
    listroi = setdiff(unique(newcci(:))',0);
    listloc =[];
    for iroi= 1:length(listroi)
        inxs = find(newcci(:)==listroi(iroi));        
        listloc = [listloc  [listroi(iroi); min(inxs)]];
    end
    [a, order]=sort(listloc(2,:),'ascend');
    for iroi= 1:length(listroi)
        inxs = find(handles.cci(:)==listloc(1,order(iroi)));
        newcci(inxs)=iroi;
    end
    selectedROIs = setdiff(unique(newcci(:))',0);
    bvalid = true;
    save(handles.fnsave,'selectedROIs','newcci','bvalid');
    close(handles.hfig);



% guidata(hObject, handles);
function ButtonDown1(hObject, eventdata, handles)
% hObject    handle to axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles = guidata(hObject);
    selectedROIs = handles.selectedROIs;    
    newcci = handles.cci;



    
    if handles.bUndo
        disp('undo')
        removedROI = selectedROIs(end);
        selectedROIs = selectedROIs(1:end);
        newcci(find(newcci(:)==removedROI))=0;
        hImg = overlayImage(handles.baseimg,newcci, handles.hax,'Discrete');


        % display roi number
        haxis = gca;
        display_ROInum(haxis,newcci,newcci)

        handles.hImg2 = hImg;
        handles.selectedROIs = selectedROIs;        
        handles.cci = newcci;
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
    display_ROInum(haxis,handles.cci,handles.cci);
    
    
    
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
    newcci = handles.cci;
    startpoint = handles.startpoint;
    endpoint = pos';
    pos_rec=sort([startpoint endpoint],2);
    
    [X,Y]=meshgrid([pos_rec(1,1):pos_rec(1,2)],[pos_rec(2,1):pos_rec(2,2)]);
    inxs=sub2ind(size(newcci),Y(:),X(:));
    selROIs = setdiff(unique(newcci(inxs)),0);   
    selROIs = selROIs(:)';
                     
    selectedROIs = setdiff(selectedROIs, selROIs);  
    ids = setdiff([1:length(selectedROIs)],selectedROIs);
    if isempty(ids)
        ids = length(selectedROIs)+1;
        newcci(inxs) = ids;
    else
        newcci(inxs)= ids(1);
    end
    selectedROIs = [selectedROIs(:); ids(1)];    
    hImg = overlayImage(handles.baseimg,newcci, handles.hax,'Discrete');


    % display roi number
    haxis = gca;
    display_ROInum(haxis,newcci,newcci)



    handles.hImg2 = hImg;
    
    handles.cci = newcci;
    handles.selectedROIs = selectedROIs;    
    handles.bUndo = true;
    set(handles.hbUndo,'Enable','on')
    set(h,'WindowButtonMotionFcn','')
    set(h,'WindowButtonUpFcn','') 
    guidata(h,handles);    
    
end
  

    
