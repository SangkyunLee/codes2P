function mansep(hFig,template,inx_in,inxb,inxcell)

global clststruct;
% 


inxin=unique([inx_in inxb]);
[ys xs]=ind2sub(size(template),inxin);
cx=round((min(xs)+max(xs))/2);
cy=round((min(ys)+max(ys))/2);
dist=round(max(max(xs)-cx, max(ys)-cy)*1.5);
xs=[cx-dist:cx+dist];
ys=[cy-dist:cy+dist];
xs=xs(find(xs>0 & xs<=size(template,2)));
ys=ys(find(ys>0 & ys<=size(template,1)));
dsize=[length(ys),length(xs)];
a=zeros(size(template));
a(ys,xs)=1; 
inxfov=find(a==1);

b=zeros(size(template));
b(inxb)=1;
c=zeros(dsize);
c(:)=b(inxfov);
inx_auto_boundary=find(c==1);
clstimg=zeros(dsize);
clstimg(:)=template(inxfov);
hMainImg = drawImgwithboundary(hFig,clstimg,[],[],[]);
uicontrol('Style', 'pushbutton', 'String', 'CLOSE','Position', [20 20 50 40],'Callback',{@closewin,hFig}); 
uicontrol('Style', 'text', 'String', sprintf('Cell ID:%d',inxcell),'Position', [230 390 130 20],'FontSize',15); 
clststruct.sel_vox = [1:prod(size(clstimg))];
clststruct.offsetYX = [ min(ys), min(xs)];
clststruct.FOV=size(template);
clststruct.localsize = dsize;
clststruct.clstimg = clstimg;
clststruct.inx_auto_boundary = inx_auto_boundary;

prompt = {'Enter No. clusters:'};
dlg_title = 'Number of clusters';
num_lines = 1;
def = {'2'};
nofig = inputdlg(prompt,dlg_title,num_lines,def);
nofig = str2num(nofig{1});
for iclst=1:nofig
    clststruct.sel_boundvox{iclst}=[];
    clststruct.sel_invox{iclst}=[];
    clststruct.lastvox{iclst}=[];
    hFigs(iclst)=figure('Position',[100+560*(iclst) 678 560 420]); 
    clststruct.hFigs=hFigs;
    
    hImg = drawImgwithboundary(hFigs(iclst),clstimg,[],inx_auto_boundary,[]);%    
    set(hImg,'ButtonDownFcn',{@ButtonDown1,[],hFigs(iclst),clstimg,inx_auto_boundary,iclst});
     uicontrol('Style', 'pushbutton', 'String', 'Add a ROI','Position', [20 10 70 40],'Callback',{@Addaclst,iclst});   
     
    set(hFigs(iclst),'WindowButtonDownFcn',{@wbd,iclst});
end



% ---------------------------
function wbd(h,evd,iclst)
    global clststruct;
    handles =guidata(h);
    Obj = get(h);
    handles.fseltype = get(h,'SelectionType');
    hchs = get(h,'Children');
    hImg=[];
    for ih=1:length(hchs)
        if strcmp(get(hchs(ih),'Type'),'axes')
            hImg= hchs(ih);        
        end
    end
    if isempty(hImg), error('no Image handle'); end
    Xrange=get(hImg,'XLim');
    Yrange=get(hImg,'YLim');
%     disp('down')


    pos = round(get(Obj.CurrentAxes,'CurrentPoint'));
    pos = pos(1,1:2);

    if strcmp(handles.fseltype,'normal') & ~(pos(1)<Xrange(1) || pos(1)>Xrange(2) || pos(2)<Yrange(1) || pos(2)>Yrange(2))
        clstimg = clststruct.clstimg;
        dsize=size(clstimg);


        x=round(pos(1,1));
        y=round(pos(1,2));
        inx_voxel=dsize(1)*(x-1)+y;    
        clststruct.sel_boundvox{iclst} = unique([clststruct.sel_boundvox{iclst} inx_voxel]); 

        % set the new values for the WindowButtonMotionFcn and
        % WindowButtonUpFcn
        set(h,'WindowButtonMotionFcn',{@wbm,iclst})
        set(h,'WindowButtonUpFcn',{@wbu,iclst})   
        hImg = drawImgwithboundary(h,clstimg,clststruct.sel_invox{iclst},clststruct.inx_auto_boundary,clststruct.sel_boundvox{iclst});%  
        set(hImg,'ButtonDownFcn',{@ButtonDown1,[],h,clstimg,clststruct.inx_auto_boundary,iclst});
        
    end
    
    guidata(h,handles);

% ---------------------------
function wbm(h,evd,iclst)
% executes while the mouse moves

    % disp('motion')
    global clststruct;
    handles =guidata(h);
    if ~isfield(handles,'fseltype')
        handles.fseltype=[];
    end
        
    hchs = get(h,'Children');
    hImg=[];
    for ih=1:length(hchs)
        if strcmp(get(hchs(ih),'Type'),'axes')
            hImg= hchs(ih);        
        end
    end
    if isempty(hImg), error('no Image handle'); end
    Xrange=get(hImg,'XLim');
    Yrange=get(hImg,'YLim');    


    pos = get(get(h,'CurrentAxes'),'CurrentPoint');
    pos = round(pos(1,1:2));

    if strcmp(handles.fseltype,'normal') & ~(pos(1)<Xrange(1) || pos(1)>Xrange(2) || pos(2)<Yrange(1) || pos(2)>Yrange(2))
    
        clstimg = clststruct.clstimg;
        dsize=size(clstimg);
        x=round(pos(1,1));
        y=round(pos(1,2));
        inx_voxel=dsize(1)*(x-1)+y;    
        clststruct.sel_boundvox{iclst} = unique([clststruct.sel_boundvox{iclst} inx_voxel]); 
        hImg = drawImgwithboundary(h,clstimg,clststruct.sel_invox{iclst},clststruct.inx_auto_boundary,clststruct.sel_boundvox{iclst});%  
        
    end
    guidata(h,handles);


function wbu(h,evd,iclst)
% executes when the mouse button is released

%     disp('up')

    handles = guidata(h);
     global clststruct;
    handles =guidata(h);
    hchs = get(h,'Children');
    hImg=[];
    for ih=1:length(hchs)
        if strcmp(get(hchs(ih),'Type'),'axes')
            hImg= hchs(ih);        
        end
    end
    if isempty(hImg), error('no Image handle'); end
    Xrange=get(hImg,'XLim');
    Yrange=get(hImg,'YLim');    


    pos = get(get(h,'CurrentAxes'),'CurrentPoint');
    pos = round(pos(1,1:2));

    if strcmp(handles.fseltype,'normal') & ~(pos(1)<Xrange(1) || pos(1)>Xrange(2) || pos(2)<Yrange(1) || pos(2)>Yrange(2))
    
        clstimg = clststruct.clstimg;
        dsize=size(clstimg);
        x=round(pos(1,1));
        y=round(pos(1,2));
        inx_voxel=dsize(1)*(x-1)+y;    
        clststruct.sel_boundvox{iclst} = unique([clststruct.sel_boundvox{iclst} inx_voxel]); 
        hImg = drawImgwithboundary(h,clstimg,clststruct.sel_invox{iclst},clststruct.inx_auto_boundary,clststruct.sel_boundvox{iclst});%  
        set(hImg,'ButtonDownFcn',{@ButtonDown1,[],h,clstimg,clststruct.inx_auto_boundary,iclst});
        
    end
    
    
    set(h,'WindowButtonMotionFcn','')
    set(h,'WindowButtonUpFcn','') 
    guidata(h,handles)



function closewin(hObject,eventdata,hfig)
    close(hfig);


function Addaclst(hObject,eventdata, iclst)
    global clststruct;
    
    tmp1=zeros(clststruct.FOV);
    dsize=clststruct.localsize;
    tmp2 = zeros(dsize);
    tmp2(clststruct.sel_invox{iclst})=2;
    tmp2(clststruct.sel_boundvox{iclst})=1;
    offset = clststruct.offsetYX; 
    tmp1(offset(1):offset(1)+dsize(1)-1,offset(2):offset(2)+dsize(2)-1)=tmp2;
    clststruct.sel_invoxinFOV{iclst} = find(tmp1(:)==2);
    clststruct.sel_boundvoxinFOV{iclst} = find(tmp1(:)==1);
    close(clststruct.hFigs(iclst));
    disp('add cluster');


% guidata(hObject, handles);
% % --- Executes on mouse press over axes background.
 function ButtonDown1(hObject, eventdata, handles,hfig,clstimg, inx_auto_boundary,iclst)
% hObject    handle to axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global clststruct;

dsize=size(clstimg);

%--------- undo
if strcmp(get(hfig,'SelectionType'),'alt')
    if ~isempty(clststruct.sel_boundvox{iclst})
        Obj=get(hObject);
        hParent=Obj.Parent;
        pos=get(hParent,'CurrentPoint');
        x=round(pos(1,1));
        y=round(pos(1,2));
        inx_voxel=dsize(1)*(x-1)+y;
        
        clststruct.sel_boundvox{iclst} = setdiff(clststruct.sel_boundvox{iclst}, inx_voxel);
               
        hImg = drawImgwithboundary(hfig,clstimg,clststruct.sel_invox{iclst},inx_auto_boundary,clststruct.sel_boundvox{iclst});%  
        set(hImg,'ButtonDownFcn',{@ButtonDown1,guidata(hImg),hfig,clstimg,inx_auto_boundary,iclst});    
       
    end


end






 function  [hImg] = drawImgwithboundary(hFig,clstimg,inxin_vox,inx_auto_boundary,inx_manual_boundary)
 
    
    clrm=jet;    
    aa=(clstimg(:)-min(clstimg(:)))/(max(clstimg(:))-min(clstimg(:)));
    cval=round(aa*64);
    cval(find(cval==0))=1;

    cpRF=ones(prod(size(clstimg)),3);
    cpRF(:,:)=clrm(cval,:);
    cpRF(inx_auto_boundary,:)=0;
    cpRF(inx_manual_boundary,:)=1;
    if ~isempty(inxin_vox)
        cpRF(inxin_vox,:)=ones(length(inxin_vox),1)*[0.5 0.5 0.5];
    end
    cpRF=reshape(cpRF,[size(clstimg) 3]);
    figure(hFig);
%     figure(hFig,'Menubar','none');
    [dy, dx, ~] =size(cpRF);

    hImg=image(cpRF); axis equal; xlim([1 dx]); ylim([1 dy]);
   
    
    
    
    
    