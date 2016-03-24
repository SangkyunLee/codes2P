function [out opts] = cv_classification(dataset,label, opts)
% function out = cv_classification(dataset,label, opts)
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
%         - cvmode: cross-validation method ('random','kfold')        
%         - Ncv: number of cross-validation
%         - pertest: percentage of testing data in cross-validation
%         - perval: percentage of validation data in cross-validation
%         - classifier: smlr(default), LDA(lambda1s=0, lambda2s=0);
%         - bSparseW: Boolean variable for selection of a weight matrix with the highest sparseness
        
% 2014-02-03, written by Sangkyun Lee

def_opts = struct('sel_cl',[1 2],'out_ch',[],'Nch',1,...
    'lambda1s',[1e-5  1e-1 1e-0 1e1],'lambda2s',[1e-5  1e-1 1e-0 1e1],...
    'bnorm',false,'mode',1,...
    'inxsample',1,...
    'cvmode','random','Ncv',10,'pertest',0.1,'perval',0.1,...
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
perval = opts.perval;
inxsample =opts.inxsample;


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

X = X(~isnan(sum(X,2)),:);
if isempty(X)
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
    inxs_cv = get_cvinxs_rand(label,Ncv,pertest,perval);    
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
            trdat = X(~isnan(sum(X,2)),inx_tr);
            trlab = label(inx_tr);
            tedat = X(~isnan(sum(X,2)),inx_te);
            telab = label(inx_te);
            valdat=X(~isnan(sum(X,2)),inx_val); 
            vallab =label(inx_val);
            
            
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
                y0b(find(y0==cond(1)))=1;
                y0b(find(y0==cond(2)))=-1;
                y1b(find(y1==cond(1)))=1;
                y1b(find(y1==cond(2)))=-1;
                
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
                trlab1(find(trlab==evtid(1)))=1;
                trlab1(find(trlab==evtid(2)))=-1;
                
                telab1 = zeros(size(telab));                
                telab1(find(telab==evtid(1)))=1;
                telab1(find(telab==evtid(2)))=-1;                
                
                vallab1 = zeros(size(vallab));                
                vallab1(find(vallab==evtid(1)))=1;
                vallab1(find(vallab==evtid(2)))=-1;
                
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
                trlab1(find(trlab==evtid(1)))=1;
                trlab1(find(trlab==evtid(2)))=-1;
                
                telab1 = zeros(size(telab));                
                telab1(find(telab==evtid(1)))=1;
                telab1(find(telab==evtid(2)))=-1;                
                
                vallab1 = zeros(size(vallab));                
                vallab1(find(vallab==evtid(1)))=1;
                vallab1(find(vallab==evtid(2)))=-1;
                
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
                    [w,status,history] = l1_ls(dataX,labely,gam,1e-5);                                        
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
                lambda(:,cvinx) = [lambda1; lambda2]                
                
                
            %---------------
            % sparse multi- logistic regression
            % lambda optimization included
            elseif strcmp(opts.classifier,'smlr')            
                for inx_lam=1:length(lambda1s)        
                    conf.W=[];conf.btrain = true;
                    conf.lambda1=lambda1s(inx_lam);
                    conf.lambda2=lambda2s(inx_lam);
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
 
       
        
        



