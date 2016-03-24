function select_cluster(hFig,template,cci,inxboundary)


figure(hFig);
inxin=setdiff([1:prod(size(template))],inxboundary);
clrm=jet;
a=template(inxin);
aa=(a-min(a))/(max(a)-min(a));
cval=round(aa*64);
cval(find(cval==0))=1;
cpRF=zeros(prod(size(template)),3);
cpRF(inxin,:)=clrm(cval,:);
cpRF=reshape(cpRF,[size(template) 3]);
hmImg=image(cpRF); axis equal
 
set(hmImg,'ButtonDownFcn',{@ButtonDown1,[],cci});
% set(hmImg,'ButtonDownFcn',{@ButtonDown1,guidata(hmImg),cci});


% guidata(hObject, handles);
function ButtonDown1(hObject, eventdata, handles,cci)
% hObject    handle to axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    Obj=get(hObject);
    hParent=Obj.Parent;
    pos=get(hParent,'CurrentPoint');
    x=round(pos(1,1));
    y=round(pos(1,2));
    disp(['ROI XY:' num2str(x) ',  ' num2str(y), ' CCI:' num2str(cci(y,x))]);
    
