function out = get_pairproj_evt(X,events,dir,normdim)
% function out = get_pairproj_evt(X,events,dir,bnorm)
% get projection angle in given events and direction
%
% INPUT:
%     X: N(sample) x nc(cell)
%     events: Nx1 vector containing events
%     dir: projection direction
%     normdim: normalization dimension
%     
% OUTPUT:
%     out.pairproj: pair-wise correlation
%     out.inx_proj2d: index of pair in 2d projection matrix
%     out.size_proj2d: size of 2d projection matrix
%     
% 2014-02-11 Sangkyun Lee

% [~, nc] = size(X);
eventtype = unique(events);
nevt = numel(eventtype);
pairproj = cell(1,nevt);
inx4proj2d = cell(1,nevt);
size2d = cell(1,nevt);
eps = 1e-10;
for ievt = 1: nevt
    inxsample = find(events==eventtype(ievt));    
    Xevt =  X(inxsample,:);
    if nargin<4
        [y1d, inx, d] = get_pairproj(Xevt,dir);   
    else
        [y1d, inx, d] = get_pairproj(Xevt,dir,normdim);   
    end
    pairproj{ievt} =  y1d;
    inx4proj2d{ievt} = inx;
    size2d{ievt} = [d d];
end
out.pairproj = pairproj;
out.inx_proj2d = inx4proj2d;
out.size_proj2d = size2d;
    