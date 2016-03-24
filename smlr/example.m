% To run the following code, you need 
% data_de : 3D matrix [nCell(no channel) X frames X trials] 
% conditions: condition vector to contain condition of each trial

opts.sel_cl = condset(:,icomp)'; % specify condtion identities to be classified
opts.mode =1; % mean response decoder
opts.inxsample = 1;
opts.Ncv = 5;    
opts.Nch = Ncell;
opts.out_ch = [];            
opts.classifier = 'LDA';
opts.cvmode = 'kfold';            
opts.pertest = 0.2;
opts.perval = 0;
opts.lambda1s = 0;
opts.lambda2s = 0;
[out opts] = cv_classification(data_de,conditions,opts);                
deaccALL(icomp) = mean(out.acc_test); 