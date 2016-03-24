function  [inxs_cv]= get_Kfoldcvinxs(label,Ncv, bval)
% function  [inxs_cv]= get_Kfoldcvinxs(label,Ncv)
%
% Ncv-fold cross-validation
% INPUT: 
%     label: 1xN vector (labels for data)
%     Ncv: no. of cross-validation
%     bval: boolean variable for validation data
%         when bval==1, the same size of validation set as the test data
% OUTPUT:
%     inxs_cv.inxs_test
%     inxs_cv.inxs_val
%     inxs_cv.inxs_train
% written by Sangkyun Lee 2014-02-03    


% To generate avalidation set, Ncv should be greater than 2.
if bval && Ncv<3,
    error('Ncv>=3 when bval==true')
end
    

label = label(:)';

if length(label)~= Ncv
    sel_cl=unique(label);

    sel_inx={};
    len_cv=zeros(1,length(sel_cl));
    for inx_cl=1:length(sel_cl)
        sel_inx{inx_cl}=find(label==sel_cl(inx_cl));
        len_cv(inx_cl)=length(sel_inx{inx_cl});
    end
    len_part0 = floor(len_cv/Ncv);
    len_part = repmat(len_part0,[Ncv 1]);
    inx_remain = len_cv-sum(len_part);
    for inxcl=1:length(sel_cl)
        len_part(1:inx_remain(inxcl),inxcl)=len_part(1:inx_remain(inxcl),inxcl)+1;
    end


    inxs_test = cell(1,Ncv);
    inxs_val = cell(1,Ncv);
    inxs_train = cell(1,Ncv);
    
    for inx_cv = 1 : Ncv
        inxs_test{inx_cv}=[];
        inxs_val{inx_cv}=[];
        for inx_cl = 1 : length(sel_cl)
            if inx_cv==1,
                len1 = 0;
            else
                len1 = sum(len_part(1:inx_cv-1,inx_cl));
            end            
            
            subinxs_val = [];
            if inx_cv==Ncv
                subinxs_test = sel_inx{inx_cl}(len1+1 : end);    
                if bval
                    len0 = len_part(1,inx_cl);
                    subinxs_val = sel_inx{inx_cl}(1 : len0);                    
                end
            else
                len2 = sum(len_part(1:inx_cv,inx_cl));            
                len3 = sum(len_part(1:inx_cv+1,inx_cl));
            
                subinxs_test = sel_inx{inx_cl}(len1+1 : len2);    
                if bval
                    subinxs_val = sel_inx{inx_cl}(len2+1 : len3);    
                end
            end
            inxs_test{inx_cv} = [inxs_test{inx_cv} subinxs_test];  
            inxs_val{inx_cv} = [inxs_val{inx_cv} subinxs_val];  
        end
        inxs_train{inx_cv} = setdiff(1:length(label),[inxs_test{inx_cv} inxs_val{inx_cv}]);
    end
    inxs_cv.inxs_test = inxs_test;
    inxs_cv.inxs_val = inxs_val;
    inxs_cv.inxs_train = inxs_train;
else
    inxs_cv=struct([]);   
end