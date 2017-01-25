function T=get_eyeparfn(T,mp)
% function T=get_eyeparfn(T,mp)
% T=TPSession; %object from class TPsession

if nargin<2
    mp = T.session_path;
end
fn = dir(fullfile(mp,'doc','*.xls'));
ffn = fullfile(mp,'doc',fn.name);
sheet=1;
finfo = load_logfile(ffn, sheet);

eyemovdir = dir(fullfile(mp,'*eye*'));
for i = 1 : length(eyemovdir)
    if eyemovdir(i).isdir
        eyemovdir = eyemovdir(i).name;   
        break;
    end
end
if ~ischar(eyemovdir)
    error('no directory for eyemovie');
end

ns = T.no_scan;
for is  = 1 : ns
    P = T.scans(is).Params;
    F = P.files;
    F.eyemovfn = [];
    F.eyeparfn = [];
    F.eyemovdir = eyemovdir;
    snum = strsplit(F.subpath_xml,'/');
    
    ns2 = length(finfo);
    is2=0;
    for is2 = is2+1 : ns2             
        if strcmp(finfo(is2).Image_directory,snum{end})
            F.eyemovfn = finfo(is2).eye;
            fnstr = fullfile(mp,'matlab/eyepar',[F.eyemovfn '*smallfitinfo.mat']);
            fn = dir(fnstr);
            fnstr = fullfile(mp, 'matlab/eyepar', fn.name);
            if isempty(fn)
                %error(sprintf('%s: eyepar not estimated yet',fnstr))
                error('%s: eyepar not estimated yet',fnstr);
            else
                F.eyeparfn = fn.name;
            end
            break;            
        end
    end
    T.scans(is).Params.files=F;
    
end






     