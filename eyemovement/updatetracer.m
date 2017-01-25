function updatetracer(M)
%     persistent par;
%     persistent VR;
    
    if ~exist('M','var') || isempty(M), M = get(gcf, 'UserData'); end
    frame = get(M.ui.Fcursor.sliderHandle, 'Value');
    M.DISP.Cframe = frame;
    frameL = M.DISP.frameL;
    frameR = frame-frameL(1)+1; % relative frame
    
    switch M.DISP.mode
        case 'PIXEL'            
            %pixels = M.par.out.sPIX{frame};    
            pixels = M.par.fit(frame).data;
        case 'FIT'
            fit = M.par.fit(frame).par;
            dim = M.par.searchInfo.srchdim(1:2);
            if length(fit)==3
                pixels =circle(fit(1),fit(2),fit(3),dim);
            else
                pixels=[];
            end
        case 'NONE'
            pixels=[];
    end
    
    I = M.data(:,:,frameR);
    
    hax = M.ui.imgAxes;
    draw(hax,I,pixels);    
    
    x = M.DISP.frameL(:);
    y = M.DISP.PARs;
    
    hCMxy = M.ui.CMxy;        
    plot(hCMxy,repmat(x,[1 2]),y(1:2,:)'); hold(hCMxy,'on');   
    valid =y(2,frameL-frameL(1)+1)>0;
    z = y(1:2,valid);
    m1 = 0.9*min(z(:));
    ylim(hCMxy,[max(m1,0) min(1.1*max(z(:)),1000)]);
    
     
    hPsize = M.ui.Psize;        
    plot(hPsize,x,y(3,:)'); hold(hPsize,'on');   
    z = y(3,valid);
    m2 = 0.9*min(z(:));
    ylim(hPsize,[m2 min(1.1*max(z(:)),100)]);
    
    inx_inval = find(y(1,:)==0);
    if ~isempty(inx_inval)
        m1s = m1*ones(1,length(frameR));
        m2s = m2*ones(1,length(frameR));
        
        plot(hCMxy,x(inx_inval),m1s,'Marker','*','MarkerSize',10,'Color','r');    
        plot(hPsize,x(inx_inval),m2s,'Marker','*','MarkerSize',10,'Color','r');
    end
    if valid(frameR)  
        plot(hCMxy,x(frameR),y(1,frameR),'Marker','.','MarkerSize',10,'Color','k');   
        plot(hCMxy,x(frameR),y(2,frameR),'Marker','.','MarkerSize',10,'Color','k');
        plot(hPsize,x(frameR),y(3,frameR)','Marker','.','MarkerSize',10,'Color','k');
    else
        plot(hCMxy,x(frameR),m1,'Marker','.','MarkerSize',10,'Color','k');
        plot(hPsize,x(frameR),m2,'Marker','.','MarkerSize',10,'Color','k');
    end

   

    
    str = sprintf('x=%d,y=%d',round(y(1,frameR)),round(y(2,frameR)));
    set(M.ui.ha1,'String',str,'FitBoxToText','on');
    str = sprintf('r=%d',round(y(3,frameR)));    
    set(M.ui.ha2,'String',str,'FitBoxToText','on');
    
    hold(hCMxy,'off');hold(hPsize,'off');
    
    set(M.fig,'UserData',M);


end


function draw(hax,I,pixels)
    %persistent self;        
    I(pixels)=255;
    axes(hax); cla;    
    imagesc(I);
%     colormap('gray');
    
end


function pix= circle(x,y,r,dim)
% dim [row x column]
    ang = 0:0.01:2*pi; 
    xp = x + r*cos(ang);
    yp = y + r*sin(ang);
    pix = sub2ind(dim,round(yp),round(xp));
    
%     a=zeros(dim);
%     a(round(inx))=1;
%     figure; imagesc(a)
end