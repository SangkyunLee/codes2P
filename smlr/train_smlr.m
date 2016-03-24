function [W, TR_outs, VAL_outs] = train_smlr(Xtr, label, Xval, lambda1,lambda2, MAX_ITER) 
% multi-nomial logistic regression with elastic net

% INPUT:
%     Xtr: matrix [d(feature dimension) x n(no sample)] 
%     label: vector [n x 1]
%     lambda1: |.|^2
%     lambda2: |.|
%     MAX_ITER: maximum iteration for training
% OUTPUT:
%     W: weight matrix
%     outs: [loglikelihood; logposterior; Wdiff];
%     

%for |.|^2, lambda1
%for |.|, lambda2

% 2016-02-10 
% This is a new code, This function works with test_smlr.
% The old function name is sl_smlr


Y=zeros(length(unique(label)),length(label));
CL=unique(label);
for i=1:length(CL)
    Y(i, label==CL(i))=1;
end

[m,n] = size(Y);
[d_orig,~] = size(Xtr);


inxnz=setdiff(1:d_orig,sum(abs(Xtr),2)==0);
Xtr=Xtr(inxnz,:);
Xtr = [Xtr; ones(1,n)];
[d,~] = size(Xtr);
W = [0.01*(rand(m-1,d)-0.5); zeros(1,d)]';  %[d x m]




%------ train
if exist('smlr','file')==3
    [W L]=smlr(W,Xtr,Y,lambda1,lambda2,MAX_ITER);
else
    % constant over iterations
    Bkkc = -1/2*((m-1)/m)*(sum(Xtr.^2,2));     % [d x 1], Bkk=[Bkkc; Bkkc; ... Bkkc];[d(m-1) x 1]
    delta_kkc = -1*lambda2./(Bkkc-lambda1);   % [d x 1]     
    [W L]=sl_smlr_matlab(W,Xtr,Y,Bkkc,delta_kkc,lambda1,lambda2,MAX_ITER);
end

LL=L(1,:);
Lp=L(2,:);
Wdiff =L(3,:);


%     close all            
%     figure;plot(LL); hold on;
%     plot(Lp,'r');
%     legend('log-likelihod','log-posterior')
b=W'*Xtr;
P=bsxfun(@rdivide,exp(b),(sum(exp(b))));
[~, mi]=max(P);


Wtmp=W;
W=zeros(d_orig+1,m);
W([inxnz d_orig+1],:)=Wtmp;
TR_outs.P = P;
TR_outs.estL = mi;
TR_outs.LL = LL;
TR_outs.Lp = Lp;
TR_outs.Wdiff = Wdiff;

%-------- validation
if ~isempty(Xval)
    [n2] = size(Xval,2);
    Xval = [Xval; ones(1,n2)];
    b=W'*Xval;
    P=bsxfun(@rdivide,exp(b),(sum(exp(b))));
    [~, mi]=max(P);
    VAL_outs.P = P;
    VAL_outs.estL = mi;
else
    VAL_outs.P = [];
    VAL_outs.estL = [];
end
    


















