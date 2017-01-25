function validate_thresholds(baseimg,par)
% function validate_thresholds(baseimg,par)
% GUI for selection area with a combination of given thresholds
% 
% par.parname ={'Mean','Std'};
% par.op1 ={'<p','<p'} % level1 operator
% par.pmaps = 
% par.op2 ='&'%level2 operator
% par.fnsave=''; The default is 'temp.mat'
% output = thr1<parmap1 & thr2<parmpa2

% 08-12-2014, Sangkyun Lee



fmap = ones(size(baseimg));
parnames = par.parname;
baseimg = repmat(baseimg, [1 1 3]);
sz = get(0,'ScreenSize');
winwidth = 500;
winheight = 500+length(par.parname)*20;
hfig = figure('Position',[sz(3)/2-200,200,winwidth,winheight]);
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
overlayImage(baseimg,fmap, hax,'Continuous');
axis off

% retrieve figure handle
handles = guihandles(hfig);
handles.hfig = hfig;
handles.par = par;
handles.baseimg = baseimg;
handles.fmap = fmap;
handles.hax = hax;

%% add slider controls & button

uiheight = 0.06;
uigap = 0.008;
uiwidth = 0.1;
sliderwidth = 0.25;
uiheightngap = uiheight + uigap;
uiwidthngap = uiwidth + uigap;


for ipar=1:length(parnames)
    uiYpos = 1-((1-fbottom)+ipar*uiheightngap);
    text_position = round([0.1*winwidth winheight*uiYpos uiwidth*winwidth winheight*uiheight]);
    htext1 = uicontrol('Style','text','Position',text_position,'String',parnames{ipar},'FontSize',14)

    slider_position = round([(0.1+uiwidthngap)*winwidth winheight*uiYpos sliderwidth*winwidth winheight*uiheight]);
    
    Mval = max(max(par.pmaps(:,:,ipar)));
    mval = min(min(par.pmaps(:,:,ipar)));
    switch par.op1{ipar}
        case '<p'
            defv = mval;
        case '>p'
            defv = Mval;
        otherwise
            error('not defined');
    end
    
    hs = uicontrol('Style', 'slider','Min',mval,'Max',Mval,'Value',defv,'Position', slider_position);    
    editbox_position = round([slider_position(1)+(sliderwidth+uigap)*winwidth winheight*uiYpos ,winwidth*uiwidth winheight*uiheight]); 
    hedit = uicontrol('Style','edit','Position',editbox_position,'String',num2str(defv),'FontSize',12);
    
    
    handles.hs(ipar) = hs;    
    handles.thr(ipar) =defv;
    handles.hedit(ipar) = hedit;
    set(hs,'Callback', {@update_slider}); 
    set(hedit, 'Callback', {@update_value});     
end



uiYpos = 1-((1-fbottom)+(ipar+1)*uiheightngap);
button_pos1 = round([0.1*winwidth winheight*uiYpos uiwidth*winwidth winheight*uiheight]);
hbutton1=uicontrol('Style', 'pushbutton', 'String', 'EXPORT', 'Position', button_pos1); 
set(hbutton1, 'Callback', {@exportdata});     

%%


if isfield(par,'fnsave')
    handles.fnsave = par.fnsave;
end


guidata(hfig,handles);
uiwait(hfig);

% handles = guidata(hfig); 
% par = handles.par;    
% thrs = handles.thr;
% baseimg = handles.baseimg;
% fmap = handles.fmap;
% 
% THR.par = par;
% THR.thrs =thrs;
% THR.fmap =fmap;
% THR.baseimg = baseimg;


end
function exportdata(hObj, event, handles)
    handles = guidata(hObj);    
    par = handles.par;    
    thrs = handles.thr;
    baseimg = handles.baseimg;
    fmap = handles.fmap;
    
    THR.par = par;
    THR.thrs =thrs;
    THR.fmap =fmap;
    THR.baseimg = baseimg;
    
    if isfield(handles,'fnsave')
        save(handles.fnsave,'THR');
    else
        save('THRtemp.mat','THR');
    end
    assignin('base','THR',THR);


    close(handles.hfig);
    
end

    
function update_slider(hObj, event, handles)
    handles = guidata(hObj);
    sel_inxpar = find(handles.hs == hObj);    
    par =handles.par;
    
     
    sliderval = get(hObj,'Value');
    npar = length(par.parname);
    for ipar = 1:npar
        if ipar == sel_inxpar
            thr = sliderval;
        else
            thr = handles.thr(ipar);
        end
        pmap = par.pmaps(:,:,ipar);        
        inxop = strfind(par.op1{ipar},'p')-1;
        if inxop<1,
            error('operator should be ''<p or ''>p');
        end
        newpmap = eval(['thr ' par.op1{ipar}(1:inxop) ' pmap']);
        if ipar==1,
            fmap = newpmap;
        else
            fmap = eval(['fmap ' par.op2 ' newpmap']);
        end 
    end
        
        
        

    handles.thr(sel_inxpar) = sliderval;    
    str = sprintf(['%.2f'],sliderval);    
    axes(handles.hax);
    hImg = overlayImage(handles.baseimg,fmap,handles.hax,'Continuous');
    handles.fmap = fmap;  
    set(handles.hedit(sel_inxpar),'String',str,'FontSize',12);
    guidata(hObj, handles);
end


function update_value(hObj, event, handles)
    handles = guidata(hObj);
    sel_inxpar = find(handles.hedit == hObj);
    
    editval = str2num(get(handles.hedit(sel_inxpar),'String'));
    fmap = handles.fmap;

    set(handles.hedit(sel_inxpar),'String',num2str(editval));
    set(handles.hs(sel_inxpar),'Value',editval);
    par =handles.par;
    npar = length(par.parname);
    for ipar = 1:npar
        if ipar == sel_inxpar
            thr = editval;
        else
            thr = handles.thr(ipar);
        end
        pmap = par.pmaps(:,:,ipar);        
        inxop = strfind(par.op1{ipar},'p')-1;
        if inxop<1,
            error('operator should be ''<p or ''>p');
        end
        newpmap = eval(['thr ' par.op1{ipar}(1:inxop) ' pmap']);
        if ipar==1,
            fmap = newpmap;
        else
            fmap = eval(['fmap ' par.op2 ' newpmap']);
        end 
    end
        
        
        

    handles.thr(sel_inxpar) = editval;        
    axes(handles.hax);
    hImg = overlayImage(handles.baseimg,fmap,handles.hax,'Continuous');
    axis off;
    handles.fmap = fmap;      
    guidata(hObj, handles);
end


