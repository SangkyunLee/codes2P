function [out opts] = cv_classification2(dataset_train,label_train,dataset_test,label_test, opts)
% function [out opts] = cv_classification2(dataset_train,label_train,dataset_test,label_test, opts)
% This function works when the training dataset is different from the
% testing dataset
% INPUT:
%     dataset(train, test): 3D matrix [nCell(no channel) X frames X trials]
%     label(train, test): 1D vector [1 X trials]
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
%         - cvmode: cross-validation method ('random','kfold','precal')        
%         - Ncv: number of cross-validation
%         - pertest: percentage of testing data in cross-validation, 
%                    This percentage is based on # trials of testdata     
%         - perval: percentage of validation data in cross-validation
%                    This percentage is based on # trials of traindata
%         - pertrain:  percentage of training data in cross-validation
%         - classifier: smlr(default), LDA(lambda1s=0, lambda2s=0);
%         - bSparseW: Boolean variable for selection of a weight matrix with the highest sparseness
        
% 2014-02-03, written by Sangkyun Lee
% 2016-08-29, cvmode('precal') added for pre-calculated sample indexesforCV
def_opts = struct('sel_cl',[1 2],'out_ch',[],'Nch',1,...
    'lambda1s',[1e-5  1e-1 1e-0 1e1],'lambda2s',[1e-5  1e-1 1e-0 1e1],...
    'bnorm',false,'mode',1,...
    'inxsample',1,...
    'cvmode','random','Ncv',10,'pertest',0.1,'perval',0.1,'pertrain',0.8,...
    'classifier','smlr',...
    'bSparseW',true);

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
pertrain = opts.pertrain;
perval = opts.perval;
if strcmp(opts.cvmode,'kfold')
    assert(length(label_test)~=length(label_train),'# training sample should be the same as # testing sample');
    pertrain = 1-pertest-perval;    
    assert(pertrain>0 & pertrain<1,'pertrain should be (0,1)');
elseif strcmp(opts.cvmode,'random')
    % when # training samples are different from # testing samples,
    % to get the given number of training samples from get_cvinxs_rand(label_train,Ncv,pertest,perval);   
    pertest0 = 1-pertrain-perval;
    assert(pertest0>0 &pertest0<1,'pertrain+perval should be (0,1)');

end
    

inxsample =opts.inxsample;


label_train =label_train(:)';
label_test =label_test(:)';



dataset_train=dataset_train(in_ch,:,:);
dataset_test=dataset_test(in_ch,:,:);

%%
%

if opts.mode  == 1
    dataset_train = reshape(mean(dataset_train(:,inxsample,:),2),[size(dataset_train,1) size(dataset_train,3)]);
    dataset_test = reshape(mean(dataset_test(:,inxsample,:),2),[size(dataset_test,1) size(dataset_test,3)]);
elseif opts.mode == 4    
    dataset_train = reshape(dataset_train(:,inxsample,:),[size(dataset_train,1)*length(inxsample) size(dataset_train,3)]);
    dataset_test = reshape(dataset_test(:,inxsample,:),[size(dataset_test,1)*length(inxsample) size(dataset_test,3)]);
else
    error('not implemented yet');
end

%%

dataset_train = dataset_train(~isnan(sum(dataset_train,2)),:);
dataset_test = dataset_test(~isnan(sum(dataset_test,2)),:);
if isempty(dataset_train) || isempty(dataset_test)
    out.acc_train = [];
    out.acc_test = [];
    out.lambda = [];
    out.Ws = [];
    return;
end
%% classifier setting
if strcmp(opts.classifier,'LDA')
    lambda1s =0;
    lambda2s =0;
    opts.lambda1s =0;
    opts.lambda2s =0;
    perval=0; %set validation set as [];
end



%% CV mode setting
if strcmp(opts.cvmode,'random')    
    inxs_cv_train = get_cvinxs_rand(label_train,Ncv,pertest0,perval);  
    excl_sample = cell(Ncv,1);
    if max(opts.common_sampleinx_train)>max(opts.common_sampleinx_test)
        intersecmap=zeros(max(opts.common_sampleinx_train),1);
        intersecmap(opts.common_sampleinx_train)=opts.common_sampleinx_test;
    else
        intersecmap=zeros(max(opts.common_sampleinx_test),1);
        intersecmap(opts.common_sampleinx_test)=opts.common_sampleinx_train;
    end
    
    % exclude all samples to be used in training and validation sets
    Nsample_test = length(label_test);
    for icv = 1 : Ncv
        inx_sam =union(intersect(inxs_cv_train.inxs_train{icv},opts.common_sampleinx_train), ...
            intersect(inxs_cv_train.inxs_val{icv},opts.common_sampleinx_train));
        assert( length(inx_sam)< (1-pertest)*Nsample_test,'Too many samples in testing data are used for training.');
        excl_sample{icv}= intersecmap(inx_sam);
        
    end
     
    
    inxs_cv_test = get_cvinxs_rand(label_test,Ncv,pertest,perval,excl_sample);    
    
elseif strcmp(opts.cvmode,'kfold')    
    if perval>0
        bval = true;
        opts.perval = 1/Ncv;
        opts.pertest = 1/Ncv;
    else
        bval = false;
        opts.perval = 0;
        opts.pertest = 1/Ncv;
    end
    inxs_cv_train = get_Kfoldcvinxs(label_train,Ncv,bval);
    inxs_cv_test = get_Kfoldcvinxs(label_test,Ncv,bval);
    
elseif strcmp(opts.cvmode,'precal')
    inxs_cv_train = opts.inxs_cv;
    inxs_cv_test = opts.inxs_cv;
    
else
    error('Only random or kfold cross-validation mode allowed');
end
inxs_cv.inxs_train = inxs_cv_train.inxs_train;
inxs_cv.inxs_val = inxs_cv_train.inxs_val;
inxs_cv.inxs_test = inxs_cv_test.inxs_test;



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
            trdat = dataset_train(~isnan(sum(dataset_train,2)),inx_tr);
            trlab = label_train(inx_tr);
            tedat = dataset_test(~isnan(sum(dataset_test,2)),inx_te);
            telab = label_test(inx_te);
            valdat= dataset_train(~isnan(sum(dataset_train,2)),inx_val); 
            vallab =label_train(inx_val);
            
            
            if strcmp(opts.classifier,'LDA')
%                 tic
                [class,err,~,~,coeff] = classify(tedat',trdat',trlab','linear');
%                 toc
                acc_train1 = 1-err;
                acc_test1 = length(find((class'-telab)==0))/length(telab);
                W = [coeff(1,2).linear; coeff(1,2).const];
                if cvinx==1,
                    Ws = zeros([size(W) Ncv]);
                end
                Ws(:,:,cvinx)=W;
                lambda(:,cvinx) = [0; 0];
            elseif strcmp(opts.classifier,'LDA_GPU')
                X0 = trdat';
                y0 = trlab(:);
                X1 = tedat';
                y1 = telab(:);
                
                
                cond = unique(y0);
                y0b=zeros(size(y0));
                y0b(y0==cond(1))=1;
                y0b(y0==cond(2))=-1;
                y1b=zeros(size(y1));
                y1b(y1==cond(1))=1;
                y1b(y1==cond(2))=-1;
                
                X0 = gpuArray([X0 ones(size(X0,1),1)]);
                X0T = gpuArray([X0 ones(size(X0,1),1)]');
                X1 = gpuArray([X1 ones(size(X1,1),1)]);
                y0 = gpuArray(y0b(:));
                y1 = gpuArray(y1b(:));
                
                
                
                
                X02 = X0T*X0;

                Xy =X0T*y0;
                w = X02\Xy;
                y0p = sign(X0*w);

                oney0 = parallel.gpu.GPUArray.ones(1, size(y0,1));
                err_tr = oney0*(abs(y0-y0p))/2/size(y0,1);
                y1p = sign(X1*w);
                oney1 = parallel.gpu.GPUArray.ones(1, size(y1,1));
                err_te = oney1*abs(y1-y1p)/2/size(y1,1);
                
                
                acc_train1 = 1-gather(err_tr);                
                acc_test1 = 1-gather(err_te);
                W = gather(w);
                if cvinx==1,
                    Ws = zeros([size(W) Ncv]);
                end
                Ws(:,:,cvinx)=W;
                lambda(:,cvinx) = [0; 0];
                
            elseif strcmp(opts.classifier,'SL2LDA')
                evtid = unique(trlab);
                if length(evtid)~=2, error('This code allows for binary classification'); end
                trlab1 = zeros(size(trlab));                
                trlab1(trlab==evtid(1))=1;
                trlab1(trlab==evtid(2))=-1;
                
                telab1 = zeros(size(telab));                
                telab1(telab==evtid(1))=1;
                telab1(telab==evtid(2))=-1;                
                
                vallab1 = zeros(size(vallab));                
                vallab1(vallab==evtid(1))=1;
                vallab1(vallab==evtid(2))=-1;
                
                trdat1 = [trdat; ones(1,size(trdat,2))];
                M = trdat1*trdat1';
                N = trdat1*trlab1';
                tedat1 = [tedat; ones(1,size(tedat,2))];
                valdat1 = [valdat; ones(1,size(valdat,2))];
                
                lambda1s = unique(lambda1s);
                for inx_lam=1:length(lambda1s)        
                   
                    lambda1=lambda1s(inx_lam);
                    M1 = M+lambda1*eye(size(M));
                    weight = M1\N;
                    acc_train(inx_lam) = length(find((sign(weight'*trdat1)- trlab1)==0))/length(trlab1);
                    if isempty(vallab1)
                        acc_val=0.5;
                    else
                        acc_val(inx_lam) = length(find((sign(weight'*valdat1)- vallab1)==0))/length(vallab1);
                    end
                    
                end
                
                
                mi= find(acc_val == max(acc_val));            
                
                lambda1=lambda1s(mi(end));
                M1 = M + lambda1*eye(size(M));
                weight = M1\N;
                
                acc_train1 = length(find((sign(weight'*trdat1)- trlab1)==0))/length(trlab1);
                acc_test1 = length(find((sign(weight'*tedat1)- telab1)==0))/length(telab1);
                if cvinx==1,
                    Ws = zeros([size(weight) Ncv]);
                end
                Ws(:,cvinx)=weight;
                lambda(:,cvinx) = lambda1;       
                
            %elastic net style
            elseif strcmp(opts.classifier,'SL2L1LDA')                
                evtid = unique(trlab);
                if length(evtid)~=2, error('This code allows for binary classification'); end
                trlab1 = zeros(size(trlab));                
                trlab1(trlab==evtid(1))=1;
                trlab1(trlab==evtid(2))=-1;
                
                telab1 = zeros(size(telab));                
                telab1(telab==evtid(1))=1;
                telab1(telab==evtid(2))=-1;                
                
                vallab1 = zeros(size(vallab));                
                vallab1(vallab==evtid(1))=1;
                vallab1(vallab==evtid(2))=-1;
                
                trdat1 = [trdat; ones(1,size(trdat,2))];                
                tedat1 = [tedat; ones(1,size(tedat,2))];
                valdat1 = [valdat; ones(1,size(valdat,2))];
                
                weights = zeros(size(trdat1,1),length(lambda1s));
                for inx_lam=1:length(lambda1s)                   
                    lambda1=lambda1s(inx_lam);
                    lambda2=lambda2s(inx_lam);
                    
                    labely = [trlab1 zeros(1,size(trdat1,1))]';
                    dataX = [trdat1 sqrt(lambda1)*eye(size(trdat1,1))]'/sqrt(1+lambda1);
                    %dataX = [trdat1 sqrt(lambda1)*eye(size(trdat1,1))]';
                    
                    gam = lambda2/sqrt(1+lambda1);
                    %gam = lambda2;
                    w = l1_ls(dataX,labely,gam,1e-5);                                        
                    weight = sqrt(1+lambda1)*w;
                    %weight = w;
                    
                    weights(:,inx_lam) = weight;
                    acc_train(inx_lam) = length(find((sign(weight'*trdat1)- trlab1)==0))/length(trlab1);
                    acc_val(inx_lam) = length(find((sign(weight'*valdat1)- vallab1)==0))/length(vallab1);
                    
                end
                
                
                mi= find(acc_val == max(acc_val));            
                
                lambda1=lambda1s(mi(end));
                lambda2=lambda2s(mi(end));
                weight = weights(:,mi(end));
                
                acc_train1 = acc_train(mi(end));
                acc_test1 = length(find((sign(weight'*tedat1)- telab1)==0))/length(telab1);
                if cvinx==1,
                    Ws = zeros([size(weight) Ncv]);
                end
                Ws(:,cvinx)=weight;
                lambda(:,cvinx) = [lambda1; lambda2];                
                
                
            %---------------
            % sparse multi- logistic regression
            % lambda optimization included
            elseif strcmp(opts.classifier,'smlr')            
                for inx_lam=1:length(lambda1s)        
                    conf.W=[];conf.btrain = true;
                    conf.lambda1=lambda1s(inx_lam);
                    conf.lambda2=lambda2s(inx_lam);
                    conf.bpar=false;
                    [W, out_tr] = sl_smlr(trdat, trlab,conf);                             
                    
                    if length(lambda1s)>1 && ~isempty(valdat) 
                        conf.W=W;
                        conf.btrain = false;
                        [~, out_te] = sl_smlr(valdat, [],conf); 
                        acc_val1(inx_lam) = length(find((sel_cl(out_te.estL)-vallab)==0))/length(vallab);
                    else
                        acc_train1= length(find((sel_cl(out_tr.estL)-trlab)==0))/length(trlab);
                    end
                end
                
                if length(lambda1s)>1
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
                end

                conf.W=W;
                conf.btrain = false;
                [~, out_te] = sl_smlr(tedat, [],conf); 
                acc_test1 = length(find((sel_cl(out_te.estL)-telab)==0))/length(telab);
                if cvinx==1,
                    Ws = zeros([size(W) Ncv]);
                end
                Ws(:,:,cvinx)=W;
                lambda(:,cvinx) = [conf.lambda1; conf.lambda2];
            else
                error('not specified classifier');
            end
        
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
out.inxs_cv = inxs_cv;
       
        
        



