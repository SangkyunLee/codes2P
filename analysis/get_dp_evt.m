function dp = get_dp_evt(X,events, comparevts)
% function dp = get_dp_evt(X,events, comparevts)
% 
%     INPUT:
%         X: T(samples) x m(variables)
%         events: vector[Tx1]; events for all samples
%         comparevts: cell cotaining groups to be compared
%     OUTPUT:
%         dp
    


if size(X,1) ~= size(events,1)
    error('dimension mismatch:size(X,1) ~= size(events,1)');
end

dp = zeros(size(X,2),length(comparevts));
for inxc=1:length(comparevts)
    evtids = comparevts{inxc};
    if size(evtids,1)~=2,
        error('only comparison between two conditions are allowed');
    end
    X1 = [];
    X2 = [];
    for inx = 1:size(comparevts{inxc},2)
        X1 = [X1; X(events==evtids(1,inx),:)];
        X2 = [X2; X(events==evtids(2,inx),:)];
    end
    mX1 = mean(X1,1);
    mX2 = mean(X2,1);
    stdX1 = std(X1,0,1);
    stdX2 = std(X2,0,1);
    dp1 = 2*(mX1-mX2)./(stdX1+stdX2);
    dp(:,inxc) = dp1(:);
end
    
 