function display_ROInum(haxis,cci,clr)

if nargin<3,
    clr='w';
end    
axes(haxis);
Nroi = length(unique(cci(:)))-1;
listroi = setdiff(unique(cci(:)),0);
for iroi= 1:Nroi
    inxs = find(cci(:)==listroi(iroi));
    [yc,xc] = ind2sub(size(cci),inxs(round(end/2)));
    text(xc,yc,num2str(listroi(iroi)),'Color',clr);
end    