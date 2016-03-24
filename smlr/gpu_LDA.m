% function [err_tr, err_te, w] = gpu_LDA(X1, X0, X0T, y0, y1)
% function [err_tr, err_te, w] = gpu_LDA(X1,X0,y0,y1)
%
% Linear Discriminant Analysis: with/without GPU
% When used with GPU, 
% This function works with data which composes of multiple conditions
% y=Xw+w0,
% X0(training), X1(testing): sample x dim+1 (augment matrix):[X ones]
% y0(training), y1(testing): sample x 1
% y0, y1 should be composed of 1 and -1
%
% Sangkyun Lee 2015-07-17





X02 = X0T*X0;

Xy =X0T*y0;
w = X02\Xy;
y0p = sign(X0*w);

oney0 = parallel.gpu.GPUArray.ones(1, size(y0,1));
err_tr = oney0*(abs(y0-y0p))/2/size(y0,1);


y1p = sign(X1*w);
oney1 = parallel.gpu.GPUArray.ones(1, size(y1,1));
err_te = oney1*abs(y1-y1p)/2/size(y1,1);
    



