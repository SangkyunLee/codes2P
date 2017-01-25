function updatethreshold(I)
% function updatethreshold(I)
% This function is called by reest_gui

    if ~exist('I','var') || isempty(I), I = get(gcf, 'UserData'); end
    THR = get(I.ui.THRSL.sliderHandle, 'Value');
    I.estpar.thr = THR;
    
    frame = I.DISP.Cframe;
    frameL = I.DISP.frameL;
    frameR = frame-frameL(1)+1; % relative frame
    
    Img = I.data(:,:,frameR);
    
    pixels = Img<THR;    
    
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