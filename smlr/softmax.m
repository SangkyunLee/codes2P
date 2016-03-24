function out = softmax(a,b)
% function out = softmax(A,B)
% out = sign(a)max(0,|a|-b)

out = sign(a)*max([0, abs(a)-b]);

