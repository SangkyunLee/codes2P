function [data, list_self] = loadTseries(finfo)
% function [data list_self] = loadTseries(finfo)
% Input: 
%   1. finfo should be a string containing the complete path and file name
%   2. finfo should be a struct containg (dir,fnpattern, fext)
%
% Output:
%    data: image matrix (3D) x, y, time
%    list_self: selected file list based on the finfo.fnpattern.

% Sangkyun Lee 08-30-2013
% 10-14-2013 debugged by Sangkyun Lee
% 02-05-2015 loading tiff file changed by Sangkyun Lee
% 11-08-2016 bug fix when called in matlab 2015


data=uint16([]);
if ischar(finfo)    
   [pathstr, name, ext] = fileparts(finfo); 
   if strcmp(ext,'.zip')
       filenames=unzip(finfo,pathstr);
       if length(filenames)>1
           error('One single file should be in the zip file.');
       end
   elseif strcmp(ext,'.tif')
       filenames{1} = fullfile(pathstr,[name ext]);
   end

   
   
    % new tiff loading code
    FileTif=filenames{1};
    InfoImage=imfinfo(FileTif);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    data=zeros(nImage,mImage,NumberImages,'uint16');
    
    v = version;
    mv = textscan(v,'%f');
    mv = mv{1}(1);

    if mv>8
        t = Tiff(FileTif,'r');
        
        for i = 1 : NumberImages
             setDirectory(t,i);
            data(:,:,i) = read(t);
        end
        close(t)
    elseif mv == 7.13
        FileID = tifflib('open',FileTif,'r');
        rps = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);

        for i=1:NumberImages
           tifflib('setDirectory',FileID,i);
           % Go through each strip of data.
           rps = min(rps,nImage);
           for r = 1:rps:nImage
              row_inds = r:min(nImage,r+rps-1);
              stripNum = tifflib('computeStrip',FileID,r);
              data(row_inds,:,i) = tifflib('readEncodedStrip',FileID,stripNum);
           end
        end
        tifflib('close',FileID);
    else
        error('not implemented yet');
    end
    
    if strcmp(ext,'.zip')
        delete(filenames{1});
    end
    list_self={};
   
elseif isstruct(finfo)
   if exist(finfo.dir,'dir')~=7,
       %msg = sprintf('Directory (%s) does not exist',finfo.dir)
       error('Directory (%s) does not exist',finfo.dir);
   else

        fnpattern = finfo.fnpattern;
        fext = finfo.fext;

        listf = dir([finfo.dir '\*.' fext]);

        ifnn=1;

        lenp=length(fnpattern);
        list_self = cell(1,length(listf));
        for ifn = 1:length(listf)    
            fn=listf(ifn).name;    
            if length(fn)>lenp+3 && strcmp(fn(end-3:end),['.' fext])                
                lenp = length(fnpattern);
                if strcmp(fn(1:lenp),fnpattern),
                    list_self{ifnn}=fn;
                    ifnn = ifnn+1;
                end
            end
        end
        Nframe = length(list_self);
        for iframe=1:Nframe    
            fn_frame=list_self{iframe};
            fullfn = fullfile(finfo.dir,fn_frame);
            img = imread(fullfn);
            data(:,:,iframe)= img;
        end

   end
else
    error('finfo is neither string nor structure');
end