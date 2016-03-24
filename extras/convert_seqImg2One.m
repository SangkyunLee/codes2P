function finfo = convert_seqImg2One(fext,deli,bremorgf, bzip,Nsegment)   
% function finfo = convert_seqImg2One(fext,deli,bremorgf, bzip,Nsegment)   
% fext: file extension e.g., tif
% deli: delimiter 
% bremorgf: boolean for removing orginal files
% bzip: compress image
% Nsegment

% Sangkyun Lee 08-30-2013
% Nsegment added by Sangkyun Lee 07-xx-2015
finfo =struct([]);

idir=1;
while true
    [fndata1,path_data1] = uigetfile(['*.' fext],'Select the image file');
    if ~fndata1 & ~path_data1,
        break;
    end
    tokens={};
    if strcmp(fndata1(end-7:end),'.ome.tif')
        remain = fndata1(1:end-8);
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
    fnpattern

    listf = dir([path_data1 filesep '*.' fext]);
    ifnn=1;
    list_self={};
    lenp=length(fnpattern);
    for ifn = 1:length(listf)    
        fn=listf(ifn).name;    
        if length(fn)>lenp+3 & strcmp(fn(end-3:end),['.' fext])                
            lenp = length(fnpattern);
            if strcmp(fn(1:lenp),fnpattern),
                list_self{ifnn}=fn;
                ifnn = ifnn+1;
            end
        end
    end
    finfo(idir).path = path_data1;
    finfo(idir).listf = list_self;
    finfo(idir).fnpattern =fnpattern(1:end-1);
    finfo(idir).outfile = fullfile(path_data1,[fnpattern(1:end-length(deli)) ])
    finfo(idir).bremorgf = bremorgf; 
    idir = idir+1;
end



for idir = 1: length(finfo)
    listf = finfo(idir).listf;
    Nframe = length(listf);
    if nargin<5
        Nsegment = Nframe;
    end
    
    fullfn_out = finfo(idir).outfile;
    for numseg = 1:99
        fullfn_out1 = [finfo(idir).outfile '-' sprintf('%02d',numseg) '.' fext];
        if exist(fullfn_out1)
            choice = questdlg('Would you like to overwrite?', 'Writing option','Yes', 'No','Yes');
            switch choice
                case 'Yes'
                    delete(fullfn_out1);
                case 'No'
                    finfo(idir).outfile = [finfo(idir).outfile '-1']                    
            end

        end
    end
    
    numseg=0;
    for ifn = 1: Nframe
        if mod(ifn,1000)==1,
            fprintf('Nframe: %d\n',ifn);
        end
        fullfn = fullfile(finfo(idir).path,listf{ifn});        
        img = imread(fullfn);  
        if mod(ifn,Nsegment)==1,            
            numseg = numseg + 1;
            fullfn_out = [finfo(idir).outfile '-' sprintf('%02d',numseg) '.' fext];            
        end
        imwrite(img,fullfn_out,'Compression','none','WriteMode','append');
    end
    if bremorgf
        for ifn=1:Nframe    
            fullfn = fullfile(finfo(idir).path,listf{ifn});  
            delete(fullfn);
        end
    end

    for numseg1 = 1 : numseg
        fullfn_out = [finfo(idir).outfile '-' sprintf('%02d',numseg1)  '.' fext];
        if exist(fullfn_out) & bzip            
            fullfn_out_zip = [fullfn_out(1:end-4) '.zip'];
            zip(fullfn_out_zip,fullfn_out);
            delete(fullfn_out);
        end
    end
    msg = sprintf(['Writing ' finfo(idir).fnpattern ' completed!']);
    disp(msg)
        
end
