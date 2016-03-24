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
hImg = overlayImage(handles.baseimg,fmap,handles.hax,handles.mode)
% handles.fmap = fmap;
handles.thr = sliderval;


set(handles.htext,'String',str,'FontSize',12);
set(handles.hedit,'String',str,'FontSize',12);
guidata(hObj, handles);