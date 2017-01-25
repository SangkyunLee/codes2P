function [out, opts] = par_cv_classification(dataset,label, opts)
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
%         - bshuffle: shuffle trials in each dimension independently    
% 2016-02-10, written by Sangkyun Lee
% 2016-03-06, shuffling option added

def_opts = struct('sel_cl',[1 2],'out_ch',[],'Nch',1,...
    'lambda1s',[1e-5  1e-1 1e-0 1e1],'lambda2s',[1e-5  1e-1 1e-0 1e1],...
    'mode',1,...
    'inxsample',1,...
    'cvmode','random','Ncv',10,'pertest',0.1,'perval',0.1,...
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
perval = opts.perval;
inxsample =opts.inxsample;
bshuffle = opts.bshuffle;
bias = opts.bias;

label =label(:)';
sel_labinx=[];
for inx=1:length(sel_cl)
    inxs_cl=find(label==sel_cl(inx));    
    sel_labinx=[sel_labinx inxs_cl];
end
label=label(sel_labinx);


if ~iscell(in_ch)
    dataset=dataset(in_ch,:,sel_labinx);
    if opts.mode  == 1
        d0 = reshape(mean(dataset(:,inxsample,:),2),[size(dataset,1) size(dataset,3)]);
    elseif opts.mode==4    
        d0 = reshape(dataset(:,inxsample,:),[size(dataset,1)*length(inxsample) size(dataset,3)]);
    end
    

else
    get_subdat = @(inx) dataset(inx,:,sel_labinx);     
    d0 = cellfun(get_subdat,in_ch,'UniformOutput',false);
    if opts.mode  == 1
        rsh = @(d) reshape(mean(d(:,inxsample,:),2),[size(d,1) size(d,3)]);
        d0 = cellfun(rsh,d0,'UniformOutput',false);        
    elseif opts.mode==4    
        rsh=@(d) reshape(d(:,inxsample,:),[size(d,1)*length(inxsample) size(d,3)]);
        d0 = cellfun(rsh,d0,'UniformOutput',false);
    end
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
elseif strcmp(opts.cvmode,'precal')
    inxs_cv = opts.inxs_cv;
    if length(inxs_cv.inxs_train)~=opts.Ncv
        error('incorrect pre-calculated CV indexes');
    end
else
    error('Only random or kfold cross-validation mode allowed');
end
    

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


inxs_test = inxs_cv.inxs_test;  
inxs_val = inxs_cv.inxs_val; 
inxs_train = inxs_cv.inxs_train;
mode = opts.mode;
classifier = opts.classifier; 

parfor cvinx = 1:Ncv
    
    % preparation for parloop
    if iscell(d0)
        X = d0{cvinx};
    else
        X = d0;
    end
    label1 = label;
    
    if bshuffle,
        X = shuffle_trials(label1,unique(label1),X')';        
    end
    
    lambda1s = l1;
    lambda1s = lambda1s(:);
    lambda2s = l2;
    lambda2s = lambda2s(:);
    
    acc_train1 = NaN;
    acc_test1 = NaN;
 
    
    sel_cl1 = sel_cl;
    %-----------------------
    
    inx_te = inxs_test{cvinx};  
    inx_val = inxs_val{cvinx}; 
    inx_tr = inxs_train{cvinx}; 
   
    
    
    switch mode
        case {1, 4}
            trdat = X(:,inx_tr);
            trlab = label1(inx_tr);
            tedat = X(:,inx_te);
            telab = label1(inx_te);
            valdat=X(:,inx_val); 
            vallab =label1(inx_val);
            
            %---------------
            % sparse multi- logistic regression
            % lambda optimization included
            if strcmp(classifier,'smlr')            
                acc_val1=zeros(1,length(lambda1s));
                for inx_lam=1:length(lambda1s)
                    lambda1=lambda1s(inx_lam);
                    lambda2=lambda2s(inx_lam);
                    
                    %[~, ~, out_val] = train_smlr(trdat, trlab, valdat,  lambda1, lambda2,1e+4);
                    
                    [~, ~, out_val] = train_smlr2(trdat, trlab, valdat,  lambda1, lambda2, bias,1e+4);
                    
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
    acc_train(cvinx) = acc_train1;
    acc_test(cvinx) =  acc_test1;         
end
      


out.acc_train = acc_train;
out.acc_test = acc_test;
out.lambda = lambda;
out.Ws = Ws;
out.inxs_cv = inxs_cv;
        
       
        
        



