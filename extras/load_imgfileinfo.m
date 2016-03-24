function finfo = load_imgfileinfo(fext,deli,bremorgf, bzip)   
% function finfo = load_data(fext,deli,bremorgf)   
% fext: file extension e.g., tif
% deli: delimiter 
% bremorgf: boolean for removing orginal file
% bzip: compress image


finfo =struct([]);
idir=1;
while true
    [fndata1,path_data1] = uigetfile(['*.' fext],'Select the image file');
    if ~fndata1 & ~path_data1,
        break;
    end
    tokens={};
    remain = fndata1(1:end-4);

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
    finfo(idir).outfile = fullfile(path_data1,[fnpattern(1:end-length(deli)) '.' fext])
    finfo(idir).bremorgf = bremorgf; 
    idir = idir+1;
end




for idir = 1: length(finfo)
    listf = finfo(idir).listf;
    Nframe = length(listf);
    fullfn_out = finfo(idir).outfile;
    if exist(fullfn_out)
        choice = questdlg('Would you like to overwrite?', 'Writing option','Yes', 'No','Yes');
        switch choice
            case 'Yes'
                delete(fullfn_out);
            case 'No'
                fullfn_out = [fullfn_out(1:end-4) '-1' fullfn_out(end-3:end)];
        end
                
    end
    for ifn = 1: Nframe
        fullfn = fullfile(finfo(idir).path,listf{ifn});        
        img = imread(fullfn);    
        imwrite(img,fullfn_out,'Compression','none','WriteMode','append')
    end
    if bremorgf
        for iframe=1:Nframe    
            fullfn = fullfile(finfo(idir).path,listf{ifn});  
            delete(fullfn);
        end
    end

    if exist(fullfn_out) & bzip
        fullfn_out_zip = [fullfn_out(1:end-4) '.zip'];
        zip(fullfn_out_zip,fullfn_out);
        delete(fullfn_out);
    end
    msg = sprintf(['Writing ' finfo(idir).fnpattern ' completed!']);
    disp(msg)
        
end
