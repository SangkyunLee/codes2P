function [inxsample, selmask1] = select_subdata(evts,sevts)
% function [inxsample, selmask1] = select_subdata(evts,sevts)
%   INPUT: 
%       evts: T (trial)x E (types of events)
%       sevts: Ex1 cell variable, each cell component contains events to
%       select

assert(iscell(sevts),'sevts should be cell variable');
assert(size(evts,2)==length(sevts),'sorts of events should be the same as the length sevts');

nevttyp = size(evts,2);
selmask = zeros(size(evts));
for ityp = 1: nevttyp
    evtsort = sevts{ityp};    
    for ievt =1 : length(evtsort)
        selmask(evts(:,ityp)==evtsort(ievt),ityp) = true;
    end
    if ityp ==1,
        selmask1 = selmask(:,ityp);
    else
        selmask1 = selmask1 &selmask(:,ityp);
    end
end
inxsample = find(selmask1);

if nargin>2
    
end

