function data = load_data(fext,deli,Nframe)   
% modified by Sangkyun Lee 2015-07-09 for loading new pairie tif file
% (.ome.tif)

[fndata1,path_data1] = uigetfile(['*.' fext],'Select the image file');
tokens={};
if strfind(fndata1,'ome.tif')
    remain=fndata1(1:end-8);
else
    remain = fndata1(1:end-4);
end

i=1;
chknumeric=[];

while ~isempty(remain)
    [token, remain] = strtok(remain,deli);
    tokens{i}=token;
    if ~isempty(str2num(token))
        chknumeric(i)=1;
    end    
    i = i +1;
end
fnpattern=[];
for ii=1:length(chknumeric)
    if ~chknumeric(ii)
        fnpattern = [fnpattern tokens{ii} deli];
    else
        break;
    end
end

listf = dir([path_data1 '\*.' fext]);

% listf = listf(3:end);
ifnn=1;
list_self={};
lenp=length(fnpattern);
for ifn = 1:length(listf)    
    fn=listf(ifn).name;    
    if length(fn)>lenp+3 && strfind(fn(end-6:end),fext)                
        lenp = length(fnpattern);
        if strcmp(fn(1:lenp),fnpattern),
            list_self{ifnn}=fn;
            ifnn = ifnn+1;
        end
    end
end
if nargin<3,
    Nframe=ifnn-1;
end
for iframe=1:Nframe    
    fn_frame=list_self{iframe};
    fullfn = fullfile(path_data1,fn_frame);
    img = double(imread(fullfn));
    
    if iframe ==1,        
        data = single([]);
    end
    data(:,:,iframe)= single(img);
end

