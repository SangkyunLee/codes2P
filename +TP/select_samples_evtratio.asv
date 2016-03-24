function inxsample = select_samples_evtratio(evts,selevt, ratio)
%function inxsample = select_samples_evtratio(evts,selevt, ratio)
% select equal or similar samples across selected events
% this function can prevent sampling bias
% evts: T (trial) x 1
% selevt: selected events
% ratio : minimum number of samples x ratio

inxconds = cell(length(selevt),1);
for inxevt = 1 : length(selevt)
    inxconds{inxevt} = find(evts==selevt(inxevt));
end
condlen = cellfun(@length,inxconds);
nsel = floor(bsxfun(@min, condlen, min(condlen)*ratio));

inxsample = zeros(sum(nsel),1);
cumn = 0;
for inxevt = 1 : length(selevt)
    inxsample(cumn+1:sum(nsel(1:inxevt)))= inxconds{inxevt}(1:nsel(inxevt));
    cumn = cumn + nsel(inxevt);
end
inxsample = sort(inxsample);    

