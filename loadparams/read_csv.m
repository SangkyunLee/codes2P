function [head,data]=read_csv(fn, nch)
% function [head,data]=read_csv(fn, nch)
% fn: fullpath of file
% nch: number of channel
% head: head{1}-- timestamps, head{3,5,..}-- channel number
% data: data(:,1)- timestamps, data(:,2..)-- data
% written by S. L. 2017-08-29

headstr = '%s';
headstr = [headstr repmat(' %s%d,',[1 nch])];
headstr = [headstr(1:end-1) '\n'];
datstr = repmat('%f, ',[1 nch+1]);
datstr = datstr(1:end-2);

fid = fopen(fn);    
head = textscan(fid,headstr);
data = textscan(fid, datstr,'CollectOutput', 1);
fclose(fid);
data = data{1};
