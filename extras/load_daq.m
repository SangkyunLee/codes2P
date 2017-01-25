function [d, h] = load_daq(fn,nch,where)
% function [d, h] = load_daq(fn,nch,where)
% fn: filename
% nch: # channel
% where: 'neurosensory','brigham', and so on
% d: data; first column is time
% h: hearder info

[~, ~, ext] = fileparts(fn);
switch where
    case {'neurosensory'}
        if strcmp(ext,'.csv')
            fid = fopen(fn);
            fl1 = cell(1,nch+1);
            fl2 = cell(1,nch+1);
            fl1{1} = '%s ';
            fl2{1} = '%f,';
            for i = 1 : nch
                fl1{i+1}='%s%d, ';                
                fl2{i+1}=' %f,';                
            end
            fl1 = cell2mat(fl1);
            fl1 = [fl1(1:end-1) '\n'];
            fl2 = cell2mat(fl2);
            fl2 = fl2(1:end-1);
            
            h = textscan(fid,fl1);        
            d = textscan(fid, fl2,'CollectOutput', 1);
            fclose(fid); 
            
        else
            d = [];
        end
    otherwise
        error('not implmented');
end
            
        
