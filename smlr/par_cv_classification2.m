function [out, opts] = par_cv_classification2(dtr,label_train,dte,label_test, opts)
%par_cv_classification2(dtr,label_train,dte,label_test, opts)
% for parallel processing
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
%         - classifier: smlr(default), 
%         - bSparseW: Boolean variable for selection of a weight matrix with the highest sparseness
%         - bshuffle: shuffle trials in each dimension independently    
% 2016-02-10, written by Sangkyun Lee
% 2016-03-06, shuffling option added
% 2016-08-29, cvmode('precal') added for pre-calculated sample indexesforCV
def_opts = struct('sel_cl',[1 2],'out_ch',[],'Nch',1,...
    'lambda1s',[1e-5  1e-1 1e-0 1e1],'lambda2s',[1e-5  1e-1 1e-0 1e1],...
    'mode',1,...
    'inxsample',1,...
    'cvmode','random','Ncv',10,'pertest',0.1,'perval',0.1,'pertrain',0.8,...
    'classifier','smlr','bias',true,...
    'bSparseW',true,'bshuffle',false);

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

if iscell(out_ch)
    gen_inch = @(y)setdiff(1:Nch,y);
    in_ch = cellfun(gen_inch,out_ch,'UniformOutput',false);
else
    in_ch=setdiff((1:Nch),out_ch);
end


lambda1s=opts.lambda1s;
lambda2s=opts.lambda2s;
% fn_prefix=opts.fn_prefix
Ncv = opts.Ncv;
bSparseW = opts.bSparseW;
pertest = opts.pertest;
pertrain = opts.pertrain;
perval = opts.perval;
inxsample =opts.inxsample;
bshuffle = opts.bshuffle;
bias = opts.bias;

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
   
[dtr,dte]= datrsh(dtr,dte,in_ch,opts.mode,inxsample);
% if ~iscell(dtr)
%     dtr = dtr(~isnan(sum(dtr,2)),:);
%     dte = dte(~isnan(sum(dte,2)),:);
%     if isempty(dtr) || isempty(dte)
%         out.acc_train = [];
%         out.acc_test = [];
%         out.lambda = [];
%         out.Ws = [];
%         return;
%     end
% end



%% shuffle data
if bshuffle,
    if iscell(dtr)
        wr_sh = @(dtr) shuffle_trials(label_train,unique(label_train),dtr')';
        dtr = cellfun(wr_sh,dtr,'UniformOutput', false);
    else
        dtr = shuffle_trials(label_train,unique(label_train),dtr')';
        
    end

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


inxs_train = inxs_cv_train.inxs_train;
inxs_val = inxs_cv_train.inxs_val;
inxs_test = inxs_cv_test.inxs_test;

inxs_cv.inxs_train = inxs_train;
inxs_cv.inxs_val = inxs_val;
inxs_cv.inxs_test = inxs_test;


mode = opts.mode;
classifier = opts.classifier; 


acc_test=zeros(1, Ncv);
acc_train=zeros(1, Ncv);
lambda = zeros(2, Ncv);

if iscell(in_ch)
    len = length(in_ch{1});
else
    len = length(in_ch);
end

if bias
    Ws = zeros([len+1 length(sel_cl) Ncv]);
else
    Ws = zeros([len length(sel_cl) Ncv]);
end
[l1,l2]=meshgrid(lambda1s, lambda2s);

parfor cvinx = 1:Ncv
    % preparation for parloop
    if iscell(dtr)
        Xtr = dtr{cvinx};
        Xte = dte{cvinx};
    else
        Xtr = dtr;
        Xte = dte;
    end
    Ltrain = label_train;
    Ltest = label_test;
    Ltrain = Ltrain(:)';
    Ltest = Ltest(:)';
    
    acc_train1 = NaN;
    acc_test1 = NaN;
    
    lambda1s = l1;
    lambda1s = lambda1s(:);
    lambda2s = l2;
    lambda2s = lambda2s(:);
    
    
    sel_cl1 = sel_cl;
    %-----------------------------
    inx_te = inxs_test{cvinx};  
    inx_val = inxs_val{cvinx}; 
    inx_tr = inxs_train{cvinx}; 

    
    
    switch mode
        case {1, 4}
            trdat = Xtr(:,inx_tr);
            trlab = Ltrain(inx_tr);
            tedat = Xte(:,inx_te);
            telab = Ltest(inx_te);
            valdat= Xtr(:,inx_val); 
            vallab =Ltrain(inx_val);
            
            
          
                
            %---------------
            % sparse multi- logistic regression
            % lambda optimization included
            if strcmp(classifier,'smlr')   
                acc_val1=zeros(1,length(lambda1s));
                
                for inx_lam=1:length(lambda1s)
                    lambda1=lambda1s(inx_lam);
                    lambda2=lambda2s(inx_lam);
                    %[~, ~, out_val] = train_smlr(trdat, trlab, valdat,  lambda1, lambda2,1e+4);
                    [~, ~, out_val] = train_smlr2(trdat, trlab, valdat,  lambda1, lambda2,bias, 1e+4);
                    acc_val1(inx_lam) = length(find((sel_cl1(out_val.estL)-vallab)==0))/length(vallab);                      
                end
                
                
                
                % find multiple instances equal to the maximum acc.
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

                lambda1=lambda1s(mi);
                lambda2=lambda2s(mi);
                %[W, out_tr] = train_smlr(trdat, trlab,[], lambda1,lambda2, 1e+4);
                [W, out_tr] = train_smlr2(trdat, trlab,[], lambda1,lambda2,bias, 1e+4);
                acc_train1= length(find((sel_cl1(out_tr.estL)-trlab)==0))/length(trlab);


                %--- decoding acc in testing data
                [estL] = test_smlr2(tedat, W,bias);               
                acc_test1 = length(find((sel_cl1(estL)-telab)==0))/length(telab);
                
                Ws(:,:,cvinx)=W;
                lambda(:,cvinx) = [lambda1; lambda2];
            else
                error('not specified classifier');
            end
        
        otherwise
            fprintf('not implemented yet');
    end
    if isnan(acc_train1)
        acc_train(cvinx) = NaN;
        acc_test(cvinx) =  NaN;
    else
        acc_train(cvinx) = acc_train1;
        acc_test(cvinx) =  acc_test1;         
    end
end
      


out.acc_train = acc_train;
out.acc_test = acc_test;
out.lambda = lambda;
out.Ws = Ws;
out.inxs_cv = inxs_cv;
end


function [dtr,dte]= datrsh(dtr,dte,in_ch,mode,inxsample)
    if iscell(in_ch)
        get_dtr = @(inx) dtr(inx,:,:);     
        get_dte = @(inx) dte(inx,:,:);             
        dtr = cellfun(get_dtr,in_ch,'UniformOutput',false);        
        dte = cellfun(get_dte,in_ch,'UniformOutput',false);
        if mode  == 1
            rsh = @(d) reshape(mean(d(:,inxsample,:),2),[size(d,1) size(d,3)]);
            dtr = cellfun(rsh,dtr,'UniformOutput',false);        
            dte = cellfun(rsh,dte,'UniformOutput',false);        
        elseif mode==4    
            rsh=@(d) reshape(d(:,inxsample,:),[size(d,1)*length(inxsample) size(d,3)]);
            dtr = cellfun(rsh,dtr,'UniformOutput',false);
            dte = cellfun(rsh,dte,'UniformOutput',false);
        end

    else
        dtr=dtr(in_ch,:,:);
        dte=dte(in_ch,:,:);

        if mode  == 1
            dtr = reshape(mean(dtr(:,inxsample,:),2),[size(dtr,1) size(dtr,3)]);
            dte = reshape(mean(dte(:,inxsample,:),2),[size(dte,1) size(dte,3)]);
        elseif mode == 4    
            dtr = reshape(dtr(:,inxsample,:),[size(dtr,1)*length(inxsample) size(dtr,3)]);
            dte = reshape(dte(:,inxsample,:),[size(dte,1)*length(inxsample) size(dte,3)]);
        else
            error('not implemented yet');
        end
    end
end
       
        
        



