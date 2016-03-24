function [hImg, hfig] = overlayImage2(img1, img2, opts)
% function [hImg, hfig] = overlayImage2(img1, img2, opts)
% overlap images with transparency
%
% img1: base 2d gray image
% img2: overlapping 2d gray image
% opts: option for alpha and 

% Sangkyun Lee 2013-10-09
if nargin<3,
    opts=struct([]);
    alphalev = 0.1;
    colorindex = 2;
    hfig= figure;
else
    if isfield(opts,'alphalev')
        alphalev = opts.alphalev;
    else
        alphalev = 0.1;
    end
    if isfield(opts,'colorindex')
        colorindex = opts.colorindex;
    else
        colorindex =2;
    end
    if isfield(opts,'hAx')
        hAx = opts.hAx;
        axes(opts.hAx); cla        
    else
        hfig = figure;
        hAx = subplot(1,1,1);
    end
end
    
if isfield(opts,'DispRange')
    img1(img1(:)<opts.DispRange(1)) = opts.DispRange(1);
    img1(img1(:)>opts.DispRange(2)) = opts.DispRange(2);
end


[M,N,D1] = size(img1); 
[M1,N1,D2] = size(img2); 
if M~=M1 ||N~=N1,  error('image dimensions are different');end

if D1>1 || D2>1, error('images should be 2d'); end

hImg1 = imagesc(img1); colormap('gray');
axis image;
hold on;

img2 = (img2 - min(img2(:)))/(max(img2(:)) -min(img2(:)));
imgover = zeros(M,N,3);
imgover(:,:,colorindex)= img2;

hImg2 = image(imgover);
hold off
set(hImg2, 'AlphaData', alphalev*sign(img2))

if exist('hAx') &&   ~isempty(hAx)    
    set(hImg1, 'Parent',hAx);
    set(hImg2, 'Parent',hAx);
end
hImg = [hImg1 hImg2];
if ~exist('hfig')
    hfig=[];
end
    




