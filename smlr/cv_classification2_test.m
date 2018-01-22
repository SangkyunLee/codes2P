function [out, opts] = cv_classification2_test(dtr,label_train,dte,label_test, WEIGHT, opts)
% function [out, opts] = cv_classification2_test(dtr,label_train,dte,label_test, WEIGHT, opts)
% This function should co-work with either par_cv_classification2.m  or cv_classification2.m 
% This function works when the training dataset is different from the
% testing dataset
% INPUT:
%     dataset(train, test): 3D matrix [nCell(no channel) X frames X trials]
%     label(train, test): 1D vector [1 X trials]
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


[dtr,dte]= datrsh(dtr,dte,in_ch,opts.mode,inxsample);




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
if strcmp(opts.cvmode,'precal')
    inxs_cv = opts.inxs_cv;
    if length(inxs_cv.inxs_train)~=opts.Ncv
        error('incorrect pre-calculated CV indexes');
    end
else
    error('Only precal cross-validation mode allowed');
end
    
inxs_train = inxs_cv.inxs_train;
inxs_test = inxs_cv.inxs_test;



mode = opts.mode;
classifier = opts.classifier; 


acc_test=zeros(1, Ncv);
acc_train=zeros(1, Ncv);



for cvinx = 1:Ncv
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
    
   
    
    sel_cl1 = sel_cl;
    %-----------------------------
    inx_te = inxs_test{cvinx};  
    inx_tr = inxs_train{cvinx}; 

    
    
    switch mode
        case {1, 4}
            trdat = Xtr(:,inx_tr);
            trlab = Ltrain(inx_tr);
            tedat = Xte(:,inx_te);
            telab = Ltest(inx_te);
                
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
       
        
        



