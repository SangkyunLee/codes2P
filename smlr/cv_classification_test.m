function [out, opts] = cv_classification_test(dataset,label, WEIGHT, opts)
% function out = cv_classification(dataset,label, opts)
% This function should co-work with either par_cv_classification.m  or cv_classification.m 
% 
% INPUT:
%     dataset: 3D matrix [nCell(no channel) X frames X trials]
%     label: 1D vector [1 X trials]
%     opts:
%         - sel_cl: selected classes
%         - out_ch: excluded channel
%         - Nch: number of channel (cell)
%         - mode: classification mode
%                 mode == 1, avg timesamples(one trial --> one prediction)
%                 mode ==2, time-dependent(frame-by-frame prediction)
%                 mode == 3, time-independent (collapse all time info)            
%                 mode == 4, spatio-temporal prediction (one trial --> one prediction)
%         - inxsample: index of samples to be used for classification.
%         - cvmode: cross-validation method ('precal')        
%         - Ncv: number of cross-validation
%         - classifier: smlr(default), 
%         - bshuffle: shuffle trials in each dimension independently    
% 2017-06-11, written by Sangkyun Lee



def_opts = struct('sel_cl',[1 2],'out_ch',[],'Nch',1,...
    'mode',1,...
    'inxsample',1,...
    'cvmode','precal','Ncv',[],...
    'classifier','smlr','bias',true,...
    'bshuffle',false);


fnms = fieldnames(def_opts);
for i=1:length(fnms),
    if ~isfield(opts,fnms{i}),
        opts.(fnms{i}) = def_opts.(fnms{i});
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


Ncv = opts.Ncv;

inxsample =opts.inxsample;
bshuffle = opts.bshuffle;
bias = opts.bias;

label =label(:)';
sel_labinx = cell(1, length(sel_cl));
for inx=1:length(sel_cl)
    inxs_cl=find(label==sel_cl(inx));    
    sel_labinx{inx} = inxs_cl(:)';
end
sel_labinx = cell2mat(sel_labinx);
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
if strcmp(opts.cvmode,'precal')
    inxs_cv = opts.inxs_cv;
    if length(inxs_cv.inxs_train)~=opts.Ncv
        error('incorrect pre-calculated CV indexes');
    end
else
    error('Only precal cross-validation mode allowed');
end
    

acc_test=zeros(1, Ncv);
acc_train=zeros(1, Ncv);



inxs_test = inxs_cv.inxs_test;  
inxs_train = inxs_cv.inxs_train;
mode = opts.mode;
classifier = opts.classifier; 

for cvinx = 1:Ncv
    
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
    

    acc_train1 = NaN;
    acc_test1 = NaN;
 
    
    sel_cl1 = sel_cl;
    %-----------------------
    
    inx_te = inxs_test{cvinx};  
    inx_tr = inxs_train{cvinx}; 
    
    switch mode
        case {1, 4}
            trdat = X(:,inx_tr);
            trlab = label1(inx_tr);
            tedat = X(:,inx_te);
            telab = label1(inx_te);

            
            %---------------
            % sparse multi- logistic regression
            % lambda optimization included
            if strcmp(classifier,'smlr')            
                    
                W = WEIGHT(:,:,cvinx);

                [estL] = test_smlr2(trdat, W,bias);               
                acc_train1 = length(find((sel_cl1(estL)-trlab)==0))/length(trlab);

                
                %--- decoding acc in testing data
                [estL] = test_smlr2(tedat, W,bias);               
                acc_test1 = length(find((sel_cl1(estL)-telab)==0))/length(telab);
                

                %Ws(:,:,cvinx)=W;

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

       
        
        



