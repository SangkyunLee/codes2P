function [out opts] = cv_smlr(dataset,label, opts)
% function out = cv_smlr(dataset,label, opts)
% 
% INPUT:
%     dataset: 3D matrix [nCell(no channel) X frames X trials]
%     label: 1D vector [1 X trials]
%     opts:
%         - sel_cl: selected classes
%         - out_ch: excluded channel
%         - Nch: number of channel (cell)
%         - lambda1s: regularization for L2 norm
%         - lambda2s: regularization for L1 norm
%         - bnorm: boolean for z normalization
%         - mode: classification mode
%                 mode == 1, avg timesamples(one trial --> one prediction)
%                 mode ==2, time-dependent(frame-by-frame prediction)
%                 mode == 3, time-independent (collapse all time info)            
%                 mode == 4, spatio-temporal prediction (one trial --> one prediction)
%         - inxsample: index of samples to be used for classification.
%         - Ncv: number of cross-validation
%         - pertest: percentage of testing data in cross-validation
%         - perval: percentage of validation data in cross-validation
%         - bSparseW: Boolean variable for selection of a weight matrix with the highest sparseness
        
% 2014-01-24, written by Sangkyun Lee

def_opts = struct('sel_cl',[1 2],'out_ch',[],'Nch',1,...
    'lambda1s',[1e-5  1e-1 1e-0 1e1],'lambda2s',[1e-5  1e-1 1e-0 1e1],...
    'bnorm',false,'mode',1,...
    'inxsample',1,'cvmode','random','Ncv',10,'pertest',0.1,'perval',0.1,'bSparseW',true);

if nargin < 3,
	opts = def_opts;
else
	fnms = fieldnames(def_opts);
	for i=1:length(fnms),
		if ~isfield(opts,fnms{i}),
			opts.(fnms{i}) = def_opts.(fnms{i});
		end;
	end;
end;

opts

sel_cl=opts.sel_cl;
out_ch=opts.out_ch;
Nch=opts.Nch;
in_ch=setdiff((1:Nch),out_ch);
lambda1s=opts.lambda1s;
lambda2s=opts.lambda2s;
% fn_prefix=opts.fn_prefix
Ncv = opts.Ncv;
bSparseW = opts.bSparseW;
pertest = opts.pertest;
perval = opts.perval;
inxsample =opts.inxsample;

% if isempty(opts.fn_data)
%     error('data file name is required');
% else
% %%data saved from save_shuffle_data.m
%     load([opts.fn_data])
% end


label =label(:)';
sel_labinx=[];
for inx=1:length(sel_cl)
    inxs_cl=find(label==sel_cl(inx));    
    sel_labinx=[sel_labinx inxs_cl];
end
label=label(sel_labinx);
dataset=dataset(in_ch,:,sel_labinx);


%%
%

if opts.mode  == 1
    dataset=reshape(mean(dataset(:,inxsample,:),2),[size(dataset,1) size(dataset,3)]);
elseif opts.mode==4    
    dataset=reshape(dataset(:,inxsample,:),[size(dataset,1)*length(inxsample) size(dataset,3)]);
end







%%
% for z-normalization
if opts.bnorm==true,
    % mode==2 or 3--> dataset should be 3dimension. Thus need to reshape
    % into 2dimension for z-normalization.
    [m n l]= size(dataset);
    if opts.mode==3,
        dataset=reshape(dataset,[size(dataset,1) size(dataset,2)*size(dataset,3)]);        
    elseif opts.mode==2,
        dataset=reshape(dataset,[size(dataset,1)*size(dataset,2) size(dataset,3)]);
    end
    X=ztrans(dataset,2);
    
    if opts.mode==2 || opts.mode==3,
        X = reshape(X,[m n l]);
    end 
else
    X = dataset;
end




%% CV mode setting
if strcmp(cvmode,'random')    
    inxs_cv = get_cvinxs_rand(label,Ncv,pertest,perval);    
elseif strcmp(cvmode,'kfold')    
    if perval>0
        bval = true;
        opts.perval = 1/Ncv;
        opts.pertest = 1/Ncv;
    else
        bval = false;
        opts.perval = 0;
        opts.pertest = 1/Ncv;
    end
    inxs_cv = get_Kfoldcvinxs(label,Ncv,bval);
else
    error('Only random or kfold cross-validation mode allowed');
end
    


acc_test=zeros(1, Ncv);
acc_train=zeros(1, Ncv);
lambda = zeros(2, Ncv);


[l1,l2]=meshgrid(lambda1s, lambda2s);
lambda1s=l1(:);
lambda2s=l2(:);
for cvinx = 1:Ncv
    
    inx_te=inxs_cv.inxs_test{cvinx};  
    inx_val=inxs_cv.inxs_val{cvinx}; 
    inx_tr=inxs_cv.inxs_train{cvinx}; 

    acc_val1=zeros(1,length(lambda1s));
    
    switch opts.mode
        case {1, 4}
            trdat = X(:,inx_tr);
            trlab = label(inx_tr);
            tedat = X(:,inx_te);
            telab = label(inx_te);
            valdat=X(:,inx_val); 
            vallab =label(inx_val);
            
            for inx_lam=1:length(lambda1s)        
                conf.W=[];conf.btrain = true;
                conf.lambda1=lambda1s(inx_lam);
                conf.lambda2=lambda2s(inx_lam);
                W = sl_smlr(trdat, trlab,conf);              
                

                conf.W=W;
                conf.btrain = false;
                [~, out_te] = sl_smlr(valdat, [],conf); 
                acc_val1(inx_lam) = length(find((sel_cl(out_te.estL)-vallab)==0))/length(vallab);
            end
            mi= find(max(acc_val1)==acc_val1);            
            
            % find the most sparse weight matrix with the highest accuracy
            if bSparseW                
                if length(mi)>1
                    inxes = find(lambda2s(mi)==max(lambda2s(mi)));
                    mi = mi(inxes(end));
                end
            else
                inxes = find(lambda2s(mi)==min(lambda2s(mi)));
                mi = mi(inxes(1));
            end

            conf.W=[];conf.btrain = true;
            conf.lambda1=lambda1s(mi);
            conf.lambda2=lambda2s(mi);
            [W, out_tr] = sl_smlr(trdat, trlab,conf);                
            acc_train1= length(find((sel_cl(out_tr.estL)-trlab)==0))/length(trlab);

            conf.W=W;
            conf.btrain = false;
            [~, out_te] = sl_smlr(tedat, [],conf); 
            acc_test1 = length(find((sel_cl(out_te.estL)-telab)==0))/length(telab);
            if cvinx==1,
                Ws = zeros([size(W) Ncv]);
            end
            Ws(:,:,cvinx)=W;
            lambda(:,cvinx) = [conf.lambda1; conf.lambda2];
        
        otherwise
            fprintf('not implemented yet');
    end
    acc_train(cvinx) = acc_train1;
    acc_test(cvinx) =  acc_test1;         
end
      


out.acc_train = acc_train;
out.acc_test = acc_test;
out.lambda = lambda;
out.Ws = Ws;
 
       
        
        



