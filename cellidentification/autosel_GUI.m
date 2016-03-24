function autosel_GUI(meanimg,fmap,baseimg,mode,fnsave,figname)


[dy dx]=size(meanimg);
winheight = 700;
winwidth = 800;
hfig=figure('Position', [400 200 winwidth winheight]);
if exist('figname')
    set(hfig,'Name',figname);
end
% left=0.05; bottom=0.45; fwidth=0.45; fheight=fwidth*dy/dx;
left=0.02;  fwidth=0.38; fheight=fwidth*dy/dx;
if fheight>0.5
    fheight=0.5;
    fwidth = fheight*dx/dy;
    
end
fbottom = 1-fheight-0.01;
haxs(1)=axes('Position',[left fbottom fwidth fheight],'Parent',hfig);   
img01=(meanimg-min(meanimg(:)))/(max(meanimg(:))-min(meanimg(:)));
img01 = repmat(img01,[1 1 3]);
hImg=image(img01);
axis equal; xlim([1 dx]); ylim([1 dy]);

% left=0.55; 

fleft=0.4;  fwidth=0.6; fheight=fwidth*dy/dx;
if fheight>0.8
    fheight=0.8;
    fwidth = fheight*dx/dy;
end
fbottom = 1-fheight-0.01;
haxs(2)=axes('Position',[fleft fbottom fwidth fheight],'Parent',hfig);   
hImg2 = overlayImage(baseimg,fmap, haxs(2),mode);
axis equal; xlim([1 dx]); ylim([1 dy]);


defv=0;
uiheight = 0.06;
uigap = 0.008;

text_position=[round(0.1*winwidth) round(winheight*(1-((1-fbottom)+uiheight+uigap))) round(0.185*winwidth) round(winheight*uiheight)];
htext1 = uicontrol('Style','text','Position',text_position,'String','Threshold','FontSize',14);

slider_position=[round(0.1*winwidth) round(winheight*(1-((1-fbottom)+2*uiheight+2*uigap))) round(0.3*winwidth) round(winheight*uiheight)];
if strcmp(mode,'Discrete')
    listval = setdiff(unique(fmap(:)),0);
    numcomp=[];
    for ii=1:length(listval)
        numcomp(ii)=length(find(fmap(:)==listval(ii)));
    end
    Mval = max(numcomp);
else
    Mval = max(fmap(:));
end
hs = uicontrol('Style', 'slider','Min',0,'Max',Mval,'Value',defv,'Position', slider_position);


handles.thr_name = ' ';
textstr = sprintf([handles.thr_name '%.2f'],0);
sliderbottom = slider_position(2);
text_position = [slider_position(1)+round((0.3+uigap)*winwidth) sliderbottom ,round(winwidth*0.13) round(winheight*uiheight)]; 
htext2 = uicontrol('Style','text','Position',text_position,...
    'String',textstr,'FontSize',12);
    
editbox_position = [slider_position(1) round(sliderbottom-(uiheight+uigap)*winheight),round(winwidth*0.13) round(winheight*uiheight)];
hedit = uicontrol('Style','edit','Position',editbox_position,'String',num2str(defv),'FontSize',12);

button_pos1 = [slider_position(1)+round(winwidth*0.14) round(sliderbottom-(uiheight+uigap)*winheight),round(winwidth*0.13) round(winheight*uiheight)]; 
hbutton1=uicontrol('Style', 'pushbutton', 'String', 'SAVE', 'Position', button_pos1); 
    
button_pos2 = [button_pos1(1)+round(winwidth*0.14) round(sliderbottom-(uiheight+uigap)*winheight),round(winwidth*0.13) round(winheight*uiheight)]; 
hbutton2=uicontrol('Style', 'pushbutton', 'String', 'CLOSE', 'Position', button_pos2); 


handles = guihandles(hfig);
handles.hfig = hfig;
handles.hs = hs;
handles.winsize=[winwidth winheight];
handles.slider_position = slider_position;
handles.fmap = fmap;
handles.baseimg = baseimg;
handles.thr =defv;
handles.hax = haxs(2);
handles.hedit = hedit;
handles.mode = mode;
% handles.mode = 'Continuous'
handles.fnsave = fnsave;

handles.htext = htext2;
guidata(hfig,handles);

set(hs,'Callback', {@update_slider}); 
set(hedit, 'Callback', {@update_value}); 
set(hbutton1, 'Callback', {@savedata}); 
set(hbutton2, 'Callback', {@closef}); 
% set(hs,'Callback', {@update_slider,guidata(hfig)}); 
% set(hedit, 'Callback', {@update_value,guidata(hfig)}); 
% set(hbutton1, 'Callback', {@savedata,guidata(hfig)}); 
guidata(hfig,handles);
uiwait(hfig);


function savedata(hObj, event, handles)
    handles = guidata(hObj);
    thr = handles.thr;
    fmap = handles.fmap;
    baseimg = handles.baseimg;
    bvalid = true;
    save(handles.fnsave,'thr','fmap','baseimg','bvalid');
    close(handles.hfig);
    
function closef(hObj, event, handles)
    handles = guidata(hObj);
    bvalid = false;
    save(handles.fnsave,'bvalid');
    close(handles.hfig);

    
function update_slider(hObj, event, handles)
handles = guidata(hObj);
fmap = handles.fmap;
if strcmp(handles.mode,'Continuous')
    sliderval = get(hObj,'Value');
    str = sprintf(['%.2f'],sliderval);
    
    fmap(find(fmap<sliderval))=0;
else
    sliderval = round(get(hObj,'Value'));
    str = sprintf([ '%d'],sliderval);
    
    
    num = length(unique(fmap))-1;
    for ic=1:num
       inxp_clust = find(fmap==ic);
       if length(inxp_clust) < sliderval
           fmap(inxp_clust)=0;
       end
    end
end


axes(handles.hax);
hImg = overlayImage(handles.baseimg,fmap,handles.hax,handles.mode);
% handles.fmap = fmap;
handles.thr = sliderval;


set(handles.htext,'String',str,'FontSize',12);
set(handles.hedit,'String',str,'FontSize',12);
guidata(hObj, handles);



function update_value(hObj, event, handles)
handles = guidata(hObj);

newvalue = str2num(get(handles.hedit,'String'));
fmap = handles.fmap;

set(handles.htext,'String',num2str(newvalue));
set(handles.hs,'Value',newvalue);


if strcmp(handles.mode,'Continuous')
    
    str = sprintf([ '%.2f'],newvalue);
    
    fmap(find(fmap<newvalue))=0;
else
    
    str = sprintf([ '%d'],newvalue);
    
    
    num = length(unique(fmap))-1;
    for ic=1:num
       inxp_clust = find(fmap==ic);
       if length(inxp_clust) < round(newvalue)
           fmap(inxp_clust)=0;
       end
    end
end

axes(handles.hax);
hImg = overlayImage(handles.baseimg,fmap,handles.hax,handles.mode);
% handles.fmap = fmap;
handles.thr = newvalue;




set(handles.htext,'String',str,'FontSize',12);
set(handles.hedit,'String',str,'FontSize',12);
guidata(hObj, handles);
