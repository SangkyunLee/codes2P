function [Y Beta DCT]=applyhighDCTfilter(X,ndct)
% function [Y Beta DCT]=applyhighDCTfilter(X,ndct)
%
% INPUT
%     X: T(samples) x M(ROIs)
%     ndct: number of the maximum cycle
% 2013-11-1 Sangkyun Lee

T = size(X,1);
tc = linspace(0,2*pi,T)';    
ndct1 = 2*ndct+1;
DCT(1:T,1:ndct1) = cos(tc*(0:0.5:ndct));

Beta = DCT\X;
Y = X - DCT*Beta;

