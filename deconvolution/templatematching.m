T=3000;
sT=3000;
Z= zeros(sT,1);
y=double(F(1:T));
n=eye(sT);

C = filter(1,[1, -P.gam],n);


lambda1=1000;
lambda2 =1;
% lambda3 = lambda1*0.001;
MAX_ITER=10;
rel_tol = 0.01;
Y = double(Y);
N = size(Y,2);
a = ones(N,1);
a = a/sqrt(a'*a);

lambda1s=zeros(1,MAX_ITER);
lambda2s=zeros(1,MAX_ITER);
lamns=zeros(1,MAX_ITER);
as=zeros(N,MAX_ITER);
posts=zeros(1,MAX_ITER);
post1s=zeros(1,MAX_ITER);
ns = zeros(T,MAX_ITER);


for iter=1:MAX_ITER
    lamn = 2*lambda2/lambda1;    
    [n,status]=getn(C,a,Y,lamn,rel_tol,1);

    Cn = C*n;
    
    if N>1
        CnYYCn = sqrt(Cn'*Y*Y'*Cn);
        CnCn = Cn'*Cn;
        %lama = CnYYCn-CnCn;
        a = Y'*Cn/CnYYCn;
    end
    
    E = Y-Cn*a';    
    e2 = E(:)'*E(:);
    lambda1 = T*N/e2;
    lambda2 = T/sum(n);

    
    lamns(iter)=2*lambda2/lambda1;    
    lambda1s(iter) = lambda1;
    lambda2s(iter) = lambda2;
%     post1s(iter) = T*N/2*log(lambda1) - lambda1/2*e2;
    posts(iter) = T*N/2*log(lambda1) - lambda1/2*e2 + T*log(lambda2)- lambda2*sum(n);

    ns(:,iter)=n;
    fprintf('.');
end
fprintf('done\n');



%% deconvolution on the average of timeseries on the ROI
% T=3000;
% sT=3000;
% Z= zeros(sT,1);
% y=double(F(1:T));
% n=eye(sT);
% 
% C = filter(1,[1, -P.gam],n);
% 
% 
% lambda1=10;
% lambda2 =1;
% MAX_ITER=10;
% lambda1s=zeros(1,MAX_ITER);
% lambda2s=zeros(1,MAX_ITER);
% posts=zeros(1,MAX_ITER);
% ns = zeros(T,MAX_ITER);
% rel_tol = 0.01;
% for iter=1:MAX_ITER
% %     if iter<5,
% %         rel_tol=0.1;
% %     else
% %         rel_tol = 0.01;
% %     end
%     lambda = 2*lambda2/lambda1;
%     [x,status]=l1_ls_nonneg(C,y,lambda,rel_tol);
%     e2 = sum((y-C*x).^2);
%     lambda1 = T/e2;
%     lambda2 = T/sum(x);
%     lambda1s(iter) = lambda1;
%     lambda2s(iter) = lambda2;
%     posts(iter) = T/2*log(lambda1) - lambda1/2*e2 + T*log(lambda2)- lambda2*sum(x);
%     ns(:,iter)=x;
% end
% 
