function [mresp, err] =get_meanresp(X,events, evtids)
% function [mresp, err] =get_meanresp(X,events, evtids)
% INPUT:
%     X: Timepoints (T) x nCell (C) x trials (E)
%     events: event identity of trials
%     evtids: group of event identities (cell variable)
%     
% 2013-11-13, written by Sangkyun Lee
if size(X,3) ~= length(events)
    error('no match between data and events');
end
if iscell(evtids)
    nGrp = length(evtids);
else
    nGrp =1;
end
[T C ~]=size(X);

mresp = zeros(T,C,nGrp);
err = zeros(T,C,nGrp);
for igrp=1:nGrp
    inxtrials = [];
    for ii = 1 : length(evtids{igrp})
        inxtrials = [inxtrials; find(events(:)==evtids{igrp}(ii))];
    end
    mresp(:,:,igrp) = mean(X(:,:,inxtrials),3);
    err(:,:,igrp) = std(X(:,:,inxtrials),0,3)/sqrt(length(inxtrials));
end

    
