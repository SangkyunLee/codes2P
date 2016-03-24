function  [inxs_cv]= get_cvinxs_rand(label,Ncv,pertest,perval, excl_sample)
% function  [inxs_cv]= get_cvinxs_rand(label,Ncv,pertest,perval, excl_sample)
% generation of cvinxs with random sampling
% INPUT:
%     label
%     Ncv: no. cross-validataion
%     pertest: percent of testing data from all samples
%     perval: percent of validating data from all samples
% OUTPUT:
%     inxs_cv.inxs_test
%     inxs_cv.inxs_val
%     inxs_cv.inxs_train
% written by Sangkyun Lee 2012-xx-xx
%
% modified by Sangkyun Lee 2016-01-08
% excl_sample: samples to be excluded
%% 

if nargin>4 
   assert(iscell(excl_sample) && length(excl_sample)==Ncv,'excl_sample should be Ncv size cell variable');
else
    excl_sample = cell(Ncv,1);
end
    
    

label=label(:)';
sel_cl=unique(label);

sel_inx={};
len_cv=[];
for inx_cl=1:length(sel_cl)
    sel_inx{inx_cl}=find(label==sel_cl(inx_cl));
    len_cv(inx_cl)=length(sel_inx{inx_cl});
end
len_testpercond=round(min(len_cv)*pertest);
len_valpercond=round(min(len_cv)*perval);



for inx_cv=1:Ncv
   
    inxs_test{inx_cv}=[];
    inxs_val{inx_cv}=[];
    for inx_cl=1:length(sel_cl)
        
        
        inxs_sel = setdiff(sel_inx{inx_cl}, excl_sample{inx_cv});         
        newinxs = randperm(length(inxs_sel));
        
        
        inxs_test{inx_cv}=[inxs_test{inx_cv} inxs_sel(newinxs(1:len_testpercond))];
        inxs_val{inx_cv}=[inxs_val{inx_cv} inxs_sel(newinxs(len_testpercond+1:len_testpercond+len_valpercond))];
            
        %inxs_test{inx_cv}=[inxs_test{inx_cv} sel_inx{inx_cl}(newinxs(1:len_testpercond))];
        %inxs_val{inx_cv}=[inxs_val{inx_cv} sel_inx{inx_cl}(newinxs(len_testpercond+1:len_testpercond+len_valpercond ))];
    end    
    %A=setdiff([1:length(label)],[inxs_test{inx_cv} inxs_val{inx_cv}]);
    %B=setdiff([1:length(label)],[inxs_test{inx_cv} inxs_val{inx_cv} excl_sample{inx_cv}']);
    inxs_train{inx_cv} =setdiff([1:length(label)],[inxs_test{inx_cv} inxs_val{inx_cv} excl_sample{inx_cv}']);

end

inxs_cv.inxs_test = inxs_test;
inxs_cv.inxs_val = inxs_val;
inxs_cv.inxs_train = inxs_train;
