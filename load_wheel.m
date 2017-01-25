function [wheel, time_sec]=load_wheel(dir_DAQ,inx_wheel1)
if nargin<2
    inx_wheel1 = 3;
end
ext='csv';


f1 = dir([dir_DAQ filesep '*.' ext]);
if strcmp(ext,'csv')
    if ~isempty(f1)
        fullfn_DAQ = fullfile(dir_DAQ,f1.name);        
        fid = fopen(fullfn_DAQ);
        textscan(fid, '%s %s%d, %s%d, %s%d, %s%d\n');
        DAQ_data = textscan(fid, '%f, %f, %f, %f, %f','CollectOutput', 1);
        fclose(fid)    
        DAQ_data = DAQ_data{1};    
        time_sec = DAQ_data(:,1)/1000;   
        wheel = DAQ_data(:,inx_wheel1);
    end
end