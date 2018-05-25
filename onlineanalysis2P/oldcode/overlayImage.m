function hImg = overlayImage(img01, cci, hax,mode)
% cci: cluster image 
% img: raw image 3d image
% hax: axes handle in which image is drawn

clrm=jet;   
img02 = img01(:,:,1);

Nc = max(cci(:));
cellpixs = [];
if strcmp(mode,'Discrete')
    for ic =1 :Nc 
        inxp = find(cci(:)==ic);
        img02(inxp) = ic+1;
        cellpixs = [cellpixs; inxp(:)];    
    end
else
    inxp = find(abs(cci(:))>0);
    cellpixs = inxp(:);
    img02(inxp) = 1;
end



clustmap = img02(cellpixs);
img02(cellpixs)=0;
img02 = repmat(img02,[1 1 3]);
[ys,xs] = ind2sub(size(cci),cellpixs);
img02 = img01;
img02 = reshape(img02,[size(img02,1)*size(img02,2) size(img02,3)]);
if max(clustmap)==min(clustmap)
    a =clustmap/max(clustmap);
else
    a = (clustmap-min(clustmap))/(max(clustmap)-min(clustmap)); 
end
cval=round(a*64);
cval(find(cval==0))=1;
img02(cellpixs,:) = clrm(cval,:);
img02 = reshape(img02,[size(img01,1) size(img01,2) size(img01,3)]);
hImg=image(img02);
if ~isempty(hax)    
    set(hImg, 'Parent',hax);
end
[dy, dx, ~] =size(img02);
axis equal; xlim([1 dx]); ylim([1 dy]);