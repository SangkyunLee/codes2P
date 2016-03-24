function [out gf] = thammingfilt(data, opts)
% function out = thammingfilt(data, opts)
% temporal hamming filter
% data: N(ROIs) x T(time samples)
% opts:
%   framerate: samplingrate of data (in Hz)
%   winsize: window size(in frame,odd number for the balance) to be applied

% Sangkyun Lee 2013-10-22

framerate = opts.framerate;
winsize = opts.winsize;

hwin = floor(winsize/2);

x = (-hwin:hwin);

g2 = hamming(winsize*2+1);
gf = g2/sum(g2);
out = imfilter(double(data'), gf, 'symmetric')';

