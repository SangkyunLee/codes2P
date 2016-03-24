function [out gf] = tgaussfilt(data, opts)
% function out = tgaussfilt(data, opts)
% temporal gaussian filter
% data: N(ROIs) x T(time samples)
% opts:
%   framerate: samplingrate of data (in Hz)
%   sigma: sigma of gaussian filter (in sec)
%   winsize: window size(in frame,odd number for the balance) to be applied

% Sangkyun Lee 2013-10-09

framerate = opts.framerate;
sigma = opts.sigma*framerate;
winsize = opts.winsize;

hwin = floor(winsize/2);

x = (-hwin:hwin);
g2 =exp(-(x).^2/sigma^2);
%  g2 = hamming(winsize*2+1);
% gf = g2/sqrt(g2*g2');
gf = g2/sum(g2);
out = imfilter(double(data'), gf, 'symmetric')';

