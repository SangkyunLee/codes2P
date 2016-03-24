function finfo = load_logfile(logfn, sheet)
% function finfo = load_logfile(logfn, sheet);
% loading log info from a excel file
% Sangkyun Lee 08-30-2013

[num,txt,raw] = xlsread(logfn,sheet);

finfo=struct([]);
infos = {};
 irow=1;
for icol=1:size(raw,2)
    if ~isnan(raw{irow,icol})            
        infos{icol} = raw{irow,icol};       
    end
end

kk=1;
for irow=2:size(raw,1)    
    if  ~isnan(raw{irow,1})
        for icol=1:length(infos)
            if ~isempty(infos{icol})
                finfo(kk).(infos{icol}) = raw{irow,icol};            
            end
        end
        kk = kk+1;
    end
end