function [estL] = test_smlr(tedat, W) 
% function [TE_outs] = test_smlr(tedat, W) 
%     INPUT: tedat [d x m], d: feature dimension, m: sample dimension
%            W [(d+1) x n], n: number of conditions
%     OUTPUT: estL- class index (NOTE: This is not the direct condition)

% This is a new code, This function works with train_smlr.
% The old function name is sl_smlr

    n = size(tedat,2);
    Xte = [tedat; ones(1,n)];
    b=W'*Xte;
    P=bsxfun(@rdivide,exp(b),(sum(exp(b))));
    [~, mi]=max(P);
    estL = mi;                
           