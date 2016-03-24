function Xout = shuffletrial(X,events)
% Xout = shuffle(X,events)
% X: NxTxC (N:trial, T:time samples, C: cell)


listcond = unique(events);
Xout=[];
for icon=1: length(listcond)
    inxs_trial = find(events==listcond(icon));
    X1 =X(inxs_trial,:,:);
    [N T c]=size(X1);
    newX=[];
    for icell=1:c
        neworder=randperm(N);
        newX(1:N,:,icell)=X1(neworder,:,icell);
    end
    Xout(inxs_trial,:,:)=newX;
end


