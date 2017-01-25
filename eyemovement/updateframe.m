function updateframe(I)
% function updateframe(I)
% This function is called by reest_gui

    if ~exist('I','var') || isempty(I), I = get(gcf, 'UserData'); end
    frame = get(I.ui.FRSL.sliderHandle, 'Value');
    I.DISP.Cframe = frame;
    frameL = I.DISP.frameL;
    frameR = frame-frameL(1)+1; % relative frame
    
    Img = I.data(:,:,frameR);
    thr = I.estpar.thr;
    pixels = Img<thr;    
    
    hax = I.ui.imgAxes;
    draw(hax,Img,pixels);    
    

    
    set(I.fig,'UserData',I);
    
end

function draw(hax,I,pixels)
      
    I(pixels)=255;
    axes(hax); cla;    
    imagesc(I);
%     colormap('gray');
    
end