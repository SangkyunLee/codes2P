function reest_gui(data,frameL,hMfig)
%reest_gui(data,frameL,hMfig)
% data: VR.data(:,:,frameL);
% frameList


    %I.data=VR.data(:,:,frameL);
    %if length(frameL)>1
        %data((data(:)>150))=0;
    %end
    I.data=data;
    
        

    I.estpar.cm0 = [];
    I.estpar.thr = 10;
    I.estpar.psel=[];
    I.DISP.Cframe = frameL(1);
    I.DISP.frameL = frameL;


    Isize = size(I.data);
    figName = sprintf('Re-track pupil');
    I.fig = figure('Color', 'w', 'Name', figName, ...
                   'Units', 'pixels', 'Position', [500 500 1.1*Isize(2) 1.7*Isize(1)]);
    I.ui.imgAxes = subplot('Position', [.05 .41 .909 .5882]); 


    I.ui.THRSL = mrvSlider([.05 .24 .35 .12], 'THR', ...
        'Range', [0 255], 'IntFlag', 1, 'Value', 10, ...
        'FontSize', 11, 'Color', 'w','Callback','updatethreshold');

    I.ui.FRSL = mrvSlider([.55 .24 .35 .12], 'FRAME', ...
        'Range', [frameL(1) frameL(end)], 'IntFlag', 1, 'Value', frameL(1), ...
        'FontSize', 11, 'Color', 'w','Callback','updateframe');


    I.ui.CMtxt =  uicontrol('Style','text','Units','normalized',...
        'position',[.05 .13 .1 .05],...
        'String','CM: ','Enable','off','FontSize', 11);

    I.ui.CMedit =  uicontrol('Style','edit','Units','normalized',...
        'position',[.15 .13 .15 .05],...
        'String','','Enable','on','FontSize', 11);

    I.ui.Frtxt =  uicontrol('Style','text','Units','normalized',...
        'position',[.05 .05 .1 .05],...
        'String','FR: ','Enable','off','FontSize', 11);

    framestr = sprintf('[%d:%d]',frameL(1),frameL(end));
    I.ui.Frameedit1 =  uicontrol('Style','edit','Units','normalized',...
        'position',[.15 .05 .15 .05],...
        'String',framestr,'Enable','on','FontSize', 11);


    I.ui.REEST = uicontrol('Style','pushbutton','Units','normalized',...
        'position',[.4 .05 .2 .15],...
        'String','Re-EST','Enable','on','FontSize', 11,'Callback',{@reest,hMfig});

    I.ui.INVAL = uicontrol('Style','pushbutton','Units','normalized',...
        'position',[.65 .05 .2 .15],...
        'String','INVALID','Enable','on','FontSize', 11,'Callback',{@setinvalid,hMfig});

    %-------
    
    updatethreshold(I);
    set(I.fig,'WindowButtonDownFcn',{@wbd});
    set(I.fig,'UserData',I);
end
function wbd(h,~)
m_type = get(h, 'selectionType');

    I = get(h,'UserData');    
    [m, n,~] = size(I.data);
    % disp('down')


    pos = round(get(I.ui.imgAxes,'CurrentPoint'));
    pos = pos(1,1:2);
    
    if ~(pos(1)<1 || pos(1)>n || pos(2)<1 || pos(2)>m)
        if strcmp(m_type,'normal')
            h1 = I.ui.CMedit;
            cmstr = sprintf('%d,%d',pos(1),pos(2));
            set(h1,'String',cmstr);
            I.estpar.cm0 = pos;
            I.estpar.psel=[];
            
        elseif strcmp(m_type,'alt')
            if size(I.data,3)==1
                fmap = zeros(size(I.data));
                thr = I.estpar.thr;
                fmap(I.data(:)<thr)=1;
                Li = bwlabel(double(fmap),8);
                I.estpar.psel = [I.estpar.psel Li(pos(2),pos(1))];
            end
            
        end
    end
    
    set(h,'UserData',I);
end


function reest(h,~,hMfig)

    
    k = get(h,'Parent');
    I = get(k,'UserData');
    
    set(I.ui.REEST,'enable','off');
    set(I.ui.INVAL,'enable','off');
    
    str = get(I.ui.Frameedit1,'String');
    eval(['frames=' str]);
    estpar = I.estpar;
    Nfr = length(frames);
    %
    % frames is referrenced from M.DISP.frameL
    % framesR is referred from I.DISP.frameL
    % fr2 is the original frame index
    framesR = frames - I.DISP.frameL(1)+1; 
    
    if isempty(estpar.psel)
        opt1.thr = estpar.thr;
        opt1.winsize = 1;
        opt1.bdisp=false;
        CM0 = estpar.cm0(:)*ones(size(I.DISP.frameL));  
        out = track_rawdata2(I.data,opt1,CM0,0,framesR);   
        [fitInfo, fail] = est_puppar(I.data, out.CM,out.sPIX,150,1,framesR,3);
    else
        if size(I.data,3)>1
            return;
        end
        fmap = zeros(size(I.data));
        thr = I.estpar.thr;
        fmap(I.data(:)<thr)=1;
        Li = bwlabel(double(fmap),8);
        
        psel = estpar.psel;
        MSK = zeros(size(fmap));
        for ip = 1 : length(psel)
            inx=Li(:)==psel(ip);
            MSK(inx)=1;
        end
        
        
        
        bw=bwperim(MSK,4);
        boundary=find(bw(:));

        [y, x]= ind2sub(size(MSK), boundary);
        l0 = convhull(x,y);
        x = x(l0);
        y = y(l0);
        %figure; plot(x,y,'r.-')
        boundary = sub2ind(size(MSK),y,x); 
        try
            [z, r] = fitcircle([x y]');
            fail=false;
        catch
            fail=true;
        end

        
        fitInfo.par=[z;r];
        fitInfo.data = boundary;
        out.sPIX{framesR} = find(MSK(:)==1);       
        out.CM(:,framesR) = z;
            
    end
    
    M = get(hMfig,'USerData');
    par = M.par;
    ref0 = M.DISP.frameL(1)-1;
    fr2 = frames+ref0;
    par.out.sPIX(fr2)= out.sPIX(framesR);
    par.out.CM(:,fr2)=out.CM(:,framesR);
    for ifr = 1 : Nfr
        par.fit(fr2(ifr)).par = fitInfo(framesR(ifr)).par;
        par.fit(fr2(ifr)).data = fitInfo(framesR(ifr)).data;
        if fail(ifr)
            M.DISP.PARs(:,frames(ifr))= zeros(size(M.DISP.PARs,1),1);
        else
            M.DISP.PARs(:,frames(ifr)) = fitInfo(framesR(ifr)).par;
        end
    end
    M.par =par;    
    
    
    set(M.fig,'UserData',M);
    updatetracer(M);
    set(k,'UserData',I);
    set(I.ui.REEST,'enable','on');
    set(I.ui.INVAL,'enable','on');
    close(I.fig);
end


function setinvalid(h,~,hMfig)
    k = get(h,'Parent');
    I = get(k,'UserData');
    
    set(I.ui.REEST,'enable','off');
    set(I.ui.INVAL,'enable','off');
    str = get(I.ui.Frameedit1,'String');
    eval(['frames=' str]);
    

    Nfr = length(frames);
    M = get(hMfig,'USerData');
    par = M.par;
    ref0 = M.DISP.frameL(1)-1;
    fr2 = frames+ref0;    
    for ifr=1:Nfr
        par.fit(fr2(ifr)).par = [];
        par.fit(fr2(ifr)).data = [];
        par.out.sPIX{fr2(ifr)}=[];
        par.out.CM(:,fr2(ifr))=0;
        M.DISP.PARs(:,frames(ifr))= zeros(size(M.DISP.PARs,1),1);        
    end
    M.par =par;    
    
    
    set(M.fig,'UserData',M);
    updatetracer(M);
    set(k,'UserData',I);
    set(I.ui.REEST,'enable','on');
    set(I.ui.INVAL,'enable','on');
    close(I.fig);
end
