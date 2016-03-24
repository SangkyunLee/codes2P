function [W, outs] = sl_smlr(X, label, conf) 
% multi-nomial logistic regression with elastic net

% INPUT:
%     X: matrix [d(feature dimension) x n(no sample)] 
%     label: vector [n x 1]
%     conf.lambda1: |.|^2
%     conf.lambda2: |.|
%     conf.btrain: ture(training mode), false(testing mode)
%     conf.W: weight matrix (for testing mode, it should be 
%     conf.MAX_ITER: maximum iteration for training
% OUTPUT:
%     W: weight matrix
%     outs: [loglikelihood; logposterior; Wdiff];
%     

%for |.|^2, lambda1
%for |.|, lambda2
def_conf = struct('W',[],'lambda1',0.1,'lambda2',0.1,'MAX_ITER',1e+4,'btrain',true);
fnms = fieldnames(def_conf);
for i=1:length(fnms),
    if ~isfield(conf,fnms{i}),
        conf.(fnms{i}) = def_conf.(fnms{i});
    end;
end;

MAX_ITER = conf.MAX_ITER;
lambda1 = conf.lambda1;
lambda2 = conf.lambda2;
W = conf.W;
btrain = conf.btrain;

Y=zeros(length(unique(label)),length(label));
CL=unique(label);
for i=1:length(CL)
    inxC=find(label==CL(i));
    Y(i, inxC)=1;
end

[m,n] = size(Y);
[d_orig,n] = size(X);

if btrain
    inx0=find(sum(abs(X),2)==0);
    inxnz=setdiff([1:d_orig],inx0);
    X=X(inxnz,:);
    X = [X; ones(1,n)];
    [d,n] = size(X);
else
    X = [X; ones(1,n)];
end
% keyboard;
if btrain && isempty(W)
    W = [0.01*(rand(m-1,d)-0.5); zeros(1,d)]';  %[d x m]
elseif btrain && ~isempty(W)
    if size(W,1) ~= d_orig+1
        error('Mismatch between sizes of W and X');
    else        
        if size(W,2) ~= m, error('Mismatch between sizes of W and X'); end        
        W = W([inxnz d_orig+1],:);
    end
elseif ~btrain && isempty(W)
    error('W is required.');
end




if btrain
    if exist('smlr','file')==3
        [W L]=smlr(W,X,Y,lambda1,lambda2,MAX_ITER);
    else
        % constant over iterations
        Bkkc = -1/2*((m-1)/m)*(sum(X.^2,2));     % [d x 1], Bkk=[Bkkc; Bkkc; ... Bkkc];[d(m-1) x 1]
        delta_kkc = -1*lambda2./(Bkkc-lambda1);   % [d x 1]     
        [W L]=sl_smlr_matlab(W,X,Y,Bkkc,delta_kkc,lambda1,lambda2,MAX_ITER);
    end
    
    LL=L(1,:);
    Lp=L(2,:);
    Wdiff =L(3,:);


%     close all            
%     figure;plot(LL); hold on;
%     plot(Lp,'r');
%     legend('log-likelihod','log-posterior')
    b=W'*X;
    P=bsxfun(@rdivide,exp(b),(sum(exp(b))));
    [mv mi]=max(P);


    Wtmp=W;
    W=zeros(d_orig+1,m);
    W([inxnz d_orig+1],:)=Wtmp;
    outs.P = P;
    outs.estL = mi;
    outs.LL = LL;
    outs.Lp = Lp;
    outs.Wdiff = Wdiff;
    
else
    
    b=W'*X;
    P=bsxfun(@rdivide,exp(b),(sum(exp(b))));
    [mv mi]=max(P);
    outs.P = P;
    outs.estL = mi;
end
% imagesc(Prob_ij)

















