function [y1d, inx, d] = get_pairproj(data,dir, normdim)
% function [y1d, inx, d] = get_pairproj(data,dir, normdim)
% 
%
% INPUT: 
%     data: N(sample) x M (cell)
%     dir: direction column(1), raw(2)
%     normdim: normalization in dimension
% OUTPUT:
%     y1d: pairwise-projection coefficient
%     inx: index for 2d correlation matrix
%     d: size of 2d correlation matrix
%        
eps = 1e-10;

if nargin==3 && normdim<3
    datanorm = sqrt(sum(data.^2,normdim));
    data = bsxfun(@rdivide, data,datanorm);       
end
if dir==1,    
    y2d = data'*data;
elseif dir==2    
    y2d = data*data';
else
    error('The projection works only with dir=1 or =2');
end
d = size(y2d,1);
tmpmat1 = triu(abs(y2d)+eps,1);
inx = find(tmpmat1(:)~=0);
y1d = y2d(inx);
