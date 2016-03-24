function [Nhat, historys]= getfastspikes_dFF(dFF,sT, opts)
% function [Nhat]= getfastspikes_dFF(F,sT, opts)
% INPUT:
%     F: TxC dFF matrix
%     sT: chunck length
%     opts.dt (frame duration in sec)
%     opts.tau 1(default)

% OUTPUT:
%     Nhat
%     hisotrys

%
% 2013-11-29 Sangkyun Lee


[T nCells] = size(dFF);
Nhat=zeros(T,nCells);
nchunck = round(T/sT);

historys = cell(nchunck, nCells);

for ic=1:nCells
    Fi = dFF(:,ic);        
    
    for ichunck=1:nchunck
        tinx1 = 1+(ichunck-1)*sT;
        if ichunck>1,
            tinx1 = tinx1 -100;
        end
        if ichunck == nchunck
            tinx2 = T;
        else
            tinx2 = ichunck*sT;
        end
        Y =Fi(tinx1:tinx2);
       
        [npart,~,history] = getspikes(double(Y),opts);
        if ichunck==1,
            Nhat(tinx1:tinx2,ic) = npart(1:sT);
        
        else
            Nhat(tinx1+100:tinx2,ic) = npart(101:end);            
        end 
        historys{ichunck, ic} = history;
    end
end
