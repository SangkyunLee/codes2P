function eye = detect_slpup(eye,par)
% function detect_slpup(eye,par)
% detect/exclude small and large pupil size
% par.ls: lowbound scale, par.us: upbound scale
% IQ = y(2)-y(1);
% lb = y(1) - ls*IQ;
% ub = y(2) + us*IQ;

ls = par.ls;
us = par.us;

is = length(eye);
ra = cell(is,1);
for i = 1 : is
    ra{i} = cellfun(@mean,eye(i).r);
end

% fill outliers in invalid trials
ra = cell2mat(ra);
[lb0, ub0] = get_bound(ra,10,10);
inx_nan = find(isnan(ra));
L= length(inx_nan);
i2f = floor(L/2);
i2s = L-i2f;
k = randperm(L);
ra(inx_nan(k(1:i2f)))= ones(i2f,1)*lb0;
ra(inx_nan(k(i2f+1:end)))= ones(i2s,1)*ub0;



[lb, ub] = get_bound(ra,ls,us);


for i = 1 :is
    mr = cellfun(@mean,eye(i).r);
    iinval = ~(mr>=lb & mr<=ub);
    emptycell = cell(length(find(iinval)),1);
    eye(i).xy(iinval)= emptycell;
    eye(i).r(iinval)=emptycell;
end
end

function [lb, ub] = get_bound(X,ls,us)
    y = quantile(X,[.25 .75]); 

    IQ = y(2)-y(1);
    lb = y(1)-ls*IQ;
    ub = y(2)+us*IQ;
end
