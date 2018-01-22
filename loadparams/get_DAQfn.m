function [DAQfn, Params] = get_DAQfn(Params)
% function [DAQfn, Params] = get_DAQfn(Params)
% Params.files structure required


files = Params.files;
mainpath = files.mainpath;
if isfield(files,'DAQ_subpath')
    DAQ_subpath = files.DAQ_subpath;
else
    DAQ_subpath =[];
end

DAQ_fullpath = fullfile(mainpath,DAQ_subpath);
if ~isfield(files,'DAQ_fn') || isempty(files.DAQ_fn)    
    DAQ_fn = dir(fullfile(DAQ_fullpath,'*.csv'));
    if isempty(DAQ_fn)
        DAQ_fn =[];
    else
        DAQ_fn = DAQ_fn.name;  
        Params.files.DAQ_fn = DAQ_fn;
    end
else
    DAQ_fn = files.DAQ_fn;
end

DAQfn = fullfile(DAQ_fullpath,DAQ_fn);
if exist(DAQfn,'file')~=2,    
    error('There is no .csv file in %s\n',DAQ_fullpath);
end