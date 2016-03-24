function [out dmean dstd]= ztrans(in,direc,setN)
% [out dmean dstd]= ztrans(in,direc,setN)
% v0.01 2009-05-05 Sangkyun Lee <sangkyun.lee@student.uni-tuebingen.de>
% v0.02 2009-07-18 Sangkyun Lee <sangkyun.lee@student.uni-tuebingen.de>
% v0.03 2009-11-05 Sangkyun Lee <sangkyun.lee@student.uni-tuebingen.de>
%  add direc==4, --> just subtract the mean across time samples.
%  add direc==4, --> just subtract the mean over voxels.
if direc==1,
    dmean=mean(in,direc);
    dstd=sqrt(sum((in-ones(size(in,1),1)*dmean).^2,1)/size(in,1));
    out=(in-ones(size(in,1),1)*dmean)./(ones(size(in,1),1)*(dstd+1e-15));
elseif direc==2,
    dmean=mean(in,direc);
    dstd=sqrt(sum((in-dmean*ones(1,size(in,2))).^2,2)/size(in,2));
    out=(in-dmean*ones(1,size(in,2)))./((dstd+1e-15)*ones(1,size(in,2)));
elseif direc==3 || direc==4,
    if ~exist('setN'), error('Need image numbers for runs'); end
    if sum(setN)~=size(in,2), error('Mismatch between number of samples and sum of image numbers of runs'); end
    out=[];
    inx_start=1;
    for inx_r=1:length(setN)        
        subin=in(:,inx_start:sum(setN(1:inx_r)));

        if direc==3,
            subout=ztrans(subin,2);
        else
            subdmean=mean(subin,2);        
            subout=(subin-subdmean*ones(1,size(subin,2)));
        end
        out=[out subout];
        inx_start=inx_start+setN(inx_r);
    end
elseif direc==5,
    dmean=mean(in,1);    
    out=in-ones(size(in,1),1)*dmean;
end
