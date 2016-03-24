function M = updateROIdisp(M)
if ~exist('M','var') || isempty(M), M = get(gcf, 'UserData'); end

[ysize xsize] = size(M.data.img1);
xmove = get(M.ui.xpos.sliderHandle, 'Value');
ymove = get(M.ui.ypos.sliderHandle, 'Value');

opts.alphalev=0.1;
opts.hAx = M.ui.imgAxes;
if isfield(M.ui,'DispRange')
    opts.DispRange = M.ui.DispRange;
end
newimg = zeros(3*ysize,3*xsize);

if xmove~=0 || ymove~=0,
    Xnew = (1:xsize) + xmove + xsize-1;
    Ynew = (1:ysize) + ymove + ysize-1;
    newimg(Ynew, Xnew) = M.data.roiimg;
    newimg = newimg(ysize+1:2*ysize, xsize+1:2*xsize);
    fprintf('Xmove: %d, Ymove: %d\n',xmove, ymove);
    set(M.fig, 'UserData', M);
else
    newimg = M.data.roiimg;
end


hImg = overlayImage2(M.data.img1,newimg,opts);
M.ui.Img = hImg;
drawnow; pause(.05);
set(M.fig,'UserData',M);

nROI = length(M.data.newROI);
for iroi= 1:nROI
    if ~isempty(M.data.newROI(iroi).centerPos)
        yc = M.data.newROI(iroi).centerPos(1);    
        xc = M.data.newROI(iroi).centerPos(2);        
        if ~isfield(M.data,'btext') || M.data.btext
            text(xc,yc,M.data.newROI(iroi).name,'Color','y');    
        end
    end
end


% M.data.newroi = newimg;
