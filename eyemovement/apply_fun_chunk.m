function [out, varout] = apply_fun_chunk (hf, chunksize,out,varargin)
% function [out, varout] = apply_fun_chunk (hf, chunksize,out,varargin)
% varargin{1}:data
% varargin{2}:selframes
% varagrin{3}:CM
% varargin{4}:sPIX for est_puppar



    Nframe = size(varargin{1},3);
    nchunk = ceil(Nframe/chunksize);
    selframes = varargin{2};
    
    for ichunk = 1 : nchunk

        ixx = (ichunk-1)*chunksize + 1 : ichunk*chunksize;
        if ixx(end)>Nframe
            ixx = (ichunk-1)*chunksize + 1: Nframe;
        end    
        ixx2 = (selframes>=ixx(1) & selframes<=ixx(end));
        ixx2 = selframes(ixx2);
        ixx3 = ixx2 - ixx(1)+1;

        
        subdata = varargin{1}(:,:,ixx);
        CM1 = varargin{3}(:,ixx);



        %outtmp = track_rawdata2(subdata,opt3,InitPar,CM1,1,ixx3);   
        fnstr = func2str(hf);
        if strfind(fnstr,'track_rawdata')            
            outtmp = feval(hf,subdata,CM1,ixx3);                 
            out.sPIX(ixx2) = outtmp.sPIX(ixx3);
            out.CM(:,ixx2) = outtmp.CM(:,ixx3);
        elseif strfind(fnstr,'est_puppar')       
            if length(varargin)==4
                sPIX = varargin{4}(ixx);
            else
                sPIX=[];
            end
            [outtmp, faili] = feval(hf,subdata,CM1,sPIX,ixx3);
            out(ixx2) = outtmp(ixx3);
            varout{1}(ixx2) = faili(ixx3);
        else
        end

    end
end



