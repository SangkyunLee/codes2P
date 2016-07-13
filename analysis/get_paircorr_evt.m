function out = get_paircorr_evt(X,events,dir)
% function out = get_paircorr_evt(X,events,dir)
%
% INPUT:
%     X: N(sample) x nc(cell)
%     events: Nx1 vector containing events
%     dir: correlation direction
%     
% OUTPUT:
%     out.paircorr: pair-wise correlation
%     out.inx_corr2d: index of pair in 2d correlation matrix
%     out.size_corr2d: size of 2d corrrelation matrix
%     out.pval1d: pvalue
%     out.geomean: geometric mean firing rate
%     
% 2014-01-16 Sangkyun Lee

    
% [~, nc] = size(X);
eventtype = unique(events);
nevt = numel(eventtype);
paircorr = cell(1,nevt);
pval1d = cell(1,nevt);
inx4corr2d = cell(1,nevt);
size2d = cell(1,nevt);
% meanX = zeros(nevt,nc);
% sigcorr = zeros(nc,nc);
% noisecorr = zeros(nc,nc,nevt);
for ievt = 1: nevt
    inxsample = find(events==eventtype(ievt));    
    Xevt =  X(inxsample,:);
    [y1d, inx, d, pv1d] = get_paircorr(Xevt,dir);
     paircorr{ievt} =  y1d;
     pval1d{ievt} = pv1d;
     inx4corr2d{ievt} = inx;
     size2d{ievt} = d;
     
     mXevt = mean(Xevt,1);     
%      stdXevt = var(Xevt,0,1);     
%      mXevt = stdXevt./mXevt;
     
     geomean2d = sqrt(mXevt' * mXevt);
     geomean{ievt} = geomean2d(inx);
end
out.paircorr = paircorr;
out.inx_corr2d = inx4corr2d;
out.size_corr2d = size2d;
out.pval1d = pval1d;
% out.geomean = geomean;     