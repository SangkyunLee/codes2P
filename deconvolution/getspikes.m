function [n a history]=getspikes(Y,opts)
% function [n a history]=getspikes(Y,opts)
% 
% J(n,a) = N*T*ln(lambda1)/2 - lambda1/2*||C*n*a' -Y||^2 + T*ln(lambda2)/2 
%           - lambda2*sum(n)
%   s.t. n>=0, a'*a=1
%
% INPUT:
%      
%       Y:  Fluorescence matrix (T(time sample) x M(pixel)) 
%       opts: dt 
%             tau
%             lambda1
%             lambda2
% 
% OUTPUT:
%       n:  estimated spike rates
%       a:  spatial filter
%       history: learning history
%       
%   2013-10-24, written by Sangkyun Lee


defaults=struct('dt',0.1,'tau',1,'lambda1',1e+5,'lambda2',1,'MAX_ITER',10,'rel_tol',0.01,'ndct',0,'chunk',1000);
if nargin<2
    opts = defaults;
else
    if isstruct(opts)
        fnames = fieldnames(defaults);
        for i=1:length(fnames)
            if ~isfield(opts,fnames{i})
                opts.(fnames{i}) = defaults.(fnames{i});
            end
        end
    else
        error('opts must be a struct');
    end    
end



[T N] = size(Y);
if T > opts.chunk
    chunksize = opts.chunk;    
    offset = round(chunksize*0.05);    
    if offset<10, 
        offset=10;
    end    
    nchunk = ceil(T/chunksize);    
    sT = chunksize + 2*offset;
    if sT>T, error('chuncksize+2*offset should be less than total frames'); end
    n=eye(sT);
else
    nchunk = 1;
    chunksize = T;
    n=eye(T);    
end


gam = 1 - opts.dt/opts.tau;
C = filter(1,[1, -gam],n);


lambda1 = opts.lambda1;
lambda2 = opts.dt*opts.lambda2;
MAX_ITER = opts.MAX_ITER;
rel_tol = opts.rel_tol;

Y = double(Y);
% N = size(Y,2);
a = ones(N,1);
a = a/sqrt(a'*a);

lambda1s=zeros(1,MAX_ITER);
lambda2s=zeros(1,MAX_ITER);
lamns=zeros(1,MAX_ITER);
as=zeros(N,MAX_ITER);
posts=zeros(1,MAX_ITER);
% post1s=zeros(1,MAX_ITER);
ns = zeros(T,MAX_ITER);

bfail = false;
for iter=1:MAX_ITER
    lamn = 2*lambda2/lambda1;
    n = zeros(T,1);
    Cn = zeros(T,1);
    for ichunk = 1 : nchunk     
        tinx1 = 1+(ichunk-1)*chunksize;
        if ichunk>1,
            tinx1 = tinx1 - offset;
        end        
            
        if ichunk == nchunk
            tinx2 = T;            
            if nchunk>1 && (tinx2-tinx1)< (chunksize+offset-1)
                tinx1 = T - (chunksize+offset) +1;                
            end
        else
            tinx2 = ichunk*chunksize + offset;
        end
        
        Ysub =Y(tinx1:tinx2,:);
        
        T1 = size(Ysub,1);
        C1 = C(1:T1,1:T1);
        
        [npart, status, ~]=getn(C1,a,Ysub,lamn,rel_tol,1);
        if strcmp(status,'Failed')
            Cn = n;
            bfail = true;            
            break;
        end
        Cnpart = C1*npart;
        if nchunk == 1
            n(tinx1:tinx2) = npart;
            Cn(tinx1:tinx2,1) = Cnpart;            
        else
            if ichunk==1
                n(tinx1:tinx2-offset) = npart(1:end-offset);
                Cn(tinx1:tinx2-offset,1) = Cnpart(1:end-offset);
            elseif ichunk == nchunk
                n(tinx1+offset:end) = npart(offset+1:end);
                Cn(tinx1+offset:end,1) = Cnpart(offset+1:end);            
            else
                n(tinx1+offset:tinx2-offset) = npart(offset+1:end-offset);
                Cn(tinx1+offset:tinx2-offset,1) = Cnpart(offset+1:end-offset);                        
            end
        end
        
    end

    
    if bfail
        history.lamns = NaN;
        history.lambda1s = NaN;
        history.lambda2s = NaN;
        history.posts = NaN; 
        fprintf('Fail\n');
        return;
    else
        if N>1
            CnY = Cn'*Y;
            CnYYCn = sqrt(CnY*CnY');
            a = Y'*Cn/CnYYCn;
        end        
        
        E = Y-Cn*a';    
        e2 = E(:)'*E(:);    
        lamns(iter)=2*lambda2/lambda1;    
        lambda1s(iter) = lambda1;
        lambda2s(iter) = lambda2;
        as(:,iter) = a;
        ns(:,iter)=n;
    %     post1s(iter) = T*N/2*log(lambda1) - lambda1/2*e2;
        posts(iter) = T*N/2*log(lambda1) - lambda1/2*e2 + T*N*log(lambda2)- N*lambda2*sum(n);
    end
    
    if (iter>1 && diff(posts(iter-1:iter))<1e-3) || (iter>2 && diff(posts(iter-1:iter))<1e-3*diff(posts(1:2)))
        history.lamns = lamns(1:iter);
        history.lambda1s = lambda1s(1:iter);
        history.lambda2s = lambda2s(1:iter);
        history.posts = posts(1:iter);
%         history.ns = ns(:,1:iter);
%         history.as = as(:,1:iter);
        fprintf('.done\n');
        return;
    end
    
    %% estimate parameters
    lambda1 = T*N/e2;
%     lambda2 = T/sum(n);    

    
    
    fprintf('.');
end
history.lamns = lamns(1:iter);
history.lambda1s = lambda1s(1:iter);
history.lambda2s = lambda2s(1:iter);
history.posts = posts(1:iter);
% history.ns = ns(:,1:iter);
% history.as = as(:,1:iter);
fprintf('done\n');

