function [newCM, mp, varp]=pupc(CM,winsize)

dist = sqrt(sum(CM(1:2,:).^2,1));

nf = length(dist);
varp = zeros(nf-winsize+1,1);
mp = zeros(nf-winsize+1,1);
newCM = CM;
for ii=1: nf-winsize+1
    varp(ii) = var(dist(ii:ii+winsize-1));    
    mp(ii) = mean(dist(ii:ii+winsize-1));
    newCM(:,ii+round(winsize/2)+1) = mean(CM(:,ii:ii+winsize-1),2);
end
mp = [zeros(round(winsize/2),1); mp];
varp = [zeros(round(winsize/2),1); varp];
