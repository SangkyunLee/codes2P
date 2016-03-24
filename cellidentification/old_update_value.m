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
hImg = overlayImage(handles.baseimg,fmap,handles.hax,handles.mode)
% handles.fmap = fmap;
handles.thr = newvalue;




set(handles.htext,'String',str,'FontSize',12);
set(handles.hedit,'String',str,'FontSize',12);
guidata(hObj, handles);