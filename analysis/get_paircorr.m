function [y1d, inx, d, pv1d] = get_paircorr(data,dir)
% function [y1d, inx, d, pv1d] = get_paircorr(data,dir)
% 
% INPUT: 
%     data: N(sample) x M (cell)
%     dir: direction column(1), raw(2)
% OUTPUT:
%     y1d: pairwise-correlation
%     pv1d: p-value
%     inx: index for 2d correlation matrix
%     d: size of 2d correlation matrix
%        
eps = 1e-10;
if dir==1,
    [y2d, pv2d] = corr(data);
elseif dir==2,
    [y2d, pv2d] = corr(data');
else
    error('only direction 1 or 2 specified');
end
d = size(y2d);
tmpmat1 = triu(abs(y2d)+eps,1);
inx = find(tmpmat1(:)~=0);
y1d = y2d(inx);
pv1d = pv2d(inx);

