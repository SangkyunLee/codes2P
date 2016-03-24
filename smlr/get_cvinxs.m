function  [inxs_cv]= get_cvinxs(label,Ncv)
% function  [inxs_cv]= get_cvinxs(label,Ncv)
%
% Ncv-fold cross-validation



if length(label)~= Ncv
    sel_cl=unique(label);

    sel_inx={};
    len_cv=[];
    for inx_cl=1:length(sel_cl)
        sel_inx{inx_cl}=find(label==sel_cl(inx_cl));
        len_cv(inx_cl)=length(sel_inx{inx_cl});
    end
    len_part=round(len_cv/Ncv);


    inxs_cv={};
    for inx_cv=1:Ncv
        inxs_cv{inx_cv}=[];
        for inx_cl=1:length(sel_cl)

            len=len_part(inx_cl);

            if inx_cv==Ncv
                subinxs = sel_inx{inx_cl}((inx_cv-1)*len+1 : end);    
            else
                subinxs = sel_inx{inx_cl}((inx_cv-1)*len+1 : inx_cv*len);    
            end
            inxs_cv{inx_cv}=[inxs_cv{inx_cv} subinxs];  

        end
        if (length(inxs_cv{inx_cv}) < 0.5*sum(len_part)),
            error('re-arrange the order of cv');
        end
    end
else
    for inx_cv=1:Ncv
        inxs_cv{inx_cv}=inx_cv;  
    end
    
end