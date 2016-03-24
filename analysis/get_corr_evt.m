function out = get_corr_evt(X,events)

% INPUT:
%     X: N(sample) x nc(cell)
%     events: Nx1 vector containing events
%     
% OUTPUT:
%     out.sigcorr: signal correlation
%     out.noisecorr: noise correlation
%     
% 2013-10-29 Sangkyun Lee

    
[~, nc] = size(X);
eventtype = unique(events);
nevt = numel(eventtype);
meanX = zeros(nevt,nc);
sigcorr = zeros(nc,nc);
noisecorr = zeros(nc,nc,nevt);
for ievt = 1: nevt
    inxsample = find(events==eventtype(ievt));
    meanX (ievt,:)= mean(X(inxsample,:),1);
    
    noiseX =  X(inxsample,:);%bsxfun(@minus, X(inxsample,:), meanX(ievt,:));
    noisecorr(:,:,ievt) = corr(noiseX);
end
sigcorr(:,:) = corr(X); % signal correlation
out.sigcorr = sigcorr;
out.noisecorr = noisecorr;
        
    