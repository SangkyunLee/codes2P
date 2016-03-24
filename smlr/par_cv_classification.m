function [out opts] = par_cv_classification(dataset,label, opts)
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
    'classifier','smlr',...
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
in_ch=setdiff((1:Nch),out_ch);
lambda1s=opts.lambda1s;
lambda2s=opts.lambda2s;
% fn_prefix=opts.fn_prefix
Ncv = opts.Ncv;
bSparseW = opts.bSparseW;
pertest = opts.pertest;
perval = opts.perval;
inxsample =opts.inxsample;
bshuffle = opts.bshuffle;


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

X = dataset;
X = X(~isnan(sum(X,2)),:);

if isempty(X)
    out.acc_train = [];
    out.acc_test = [];
    out.lambda = [];
    out.Ws = [];
    return;
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


% here 
Ws = zeros([size(X,1)+1 length(sel_cl) Ncv]);


[l1,l2]=meshgrid(lambda1s, lambda2s);


inxs_test = inxs_cv.inxs_test;  
inxs_val = inxs_cv.inxs_val; 
inxs_train = inxs_cv.inxs_train;
mode = opts.mode;
classifier = opts.classifier; 

parfor cvinx = 1:Ncv
    
    % preparation for parloop
    X = dataset;
    X = X(~isnan(sum(X,2)),:);
    label1 = label;
    
    if bshuffle,
        [~, X] = shuffle_trials(label1,unique(label1),X');
        X = X';
    end
    
    lambda1s = l1;
    lambda1s = lambda1s(:);
    lambda2s = l2;
    lambda2s = lambda2s(:);
 
    
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
                    [~, ~, out_val] = train_smlr(trdat, trlab, valdat,  lambda1, lambda2,1e+4); 
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
                [W, out_tr] = train_smlr(trdat, trlab,[], lambda1,lambda2, 1e+4); 
                acc_train1= length(find((sel_cl1(out_tr.estL)-trlab)==0))/length(trlab);
                
                %--- decoding acc in testing data
                [estL] = test_smlr(tedat, W);               
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
 
       
        
        


