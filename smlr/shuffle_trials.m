function [outdata, inxtrials] = shuffle_trials(evts,evtids,data)
% function [outdata, inxtrials] = shuffle_trials(evts,evtids,data)
% This shuffle function works for shuffling with the same evtids
% INPUT: 
%        evts: [Tx1] T-trial vecotr composed of different events
%        evtids: event list to be randomly shuffled.
%        data: number of shuffled vectors / data to be shuffled
%              when data is matrix, the data size should be [TxN]
%               
% OUTPUT:
%        inxtrials: [TxN] N T-trial column vectors,which are randomly shuffled
%                   Each column vector contains original trial index after shuffling
                  
% 3/06/2016 Sangkyun Lee

uevts = unique(evts);
assert(isempty(setdiff(evtids,uevts)),'evtids should be selected on the list of evts');
T = length(evts(:));
if isscalar(data)
    N = data;    
    bmatrix = false;
else
    [T0, N]=size(data);
    assert(T0==T,'not matched data size');
    bmatrix = true;
end

outdata = zeros(size(data));
inxtrials = zeros(T,N);
for ic = 1 : N
    for ievt = evtids(:)'
        inxt = find(evts==ievt);        
        inxtrials(inxt,ic)=inxt(randperm(length(inxt)));        
    end
    if bmatrix
        outdata(:,ic)=data(inxtrials(:,ic),ic);
    end
end