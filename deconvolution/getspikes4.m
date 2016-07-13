function [n Yn a b history]=getspikes4(Y,opts)
% function [n Yn a b history]=getspikes4(Y,opts)
% 
% J(n,a) = N*T*ln(lambda1)/2 - lambda1/2*||[C*n 1]*a' -Y||^2 + T*ln(lambda2)/2 
%           - lambda2*sum(n) -lambda3*sum(a)
%   s.t. n>=0, a>=0
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
%   2014-01-21, written by Sangkyun Lee
%     I update getspike2 with Yn(dff) with spatially weight sum
%         and I rescaled n with respect to the spatial filter a
%     2015-09-08 bug fixed by Sangkyun Lee
%     When the last chunk is smaller than the side offset, error occurs.
%     I fixed it and when the last chunk is smaller than 400 frames, 
%     the last chunk is added to the (last-1) chunk.
% 2016-04-12, I change b into vector    


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
    % if the last chunk is too small (i.e.,less than 400 frames, it is
    % added in the (last -1) chunk.
    lastchunksize = mod(T,chunksize);
    if lastchunksize > 400    
        nchunk = ceil(T/chunksize);    
        sT = chunksize + 2*offset;
    else
        nchunk = floor(T/chunksize);
        sT1 = chunksize + 2* offset;
        sT2  = chunksize + offset + lastchunksize;
        if sT1>sT2
            sT = sT1;
        else
            sT =sT2;
        end
        
    end  
    
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
bs=zeros(N,MAX_ITER);
posts=zeros(1,MAX_ITER);

ns = zeros(T,MAX_ITER);
b=0;
bfail = false;
for iter=1:MAX_ITER
    lamn = 2*lambda2/lambda1;
    n = zeros(T,1);
    Cn = zeros(T,1);

    Y1 = bsxfun(@minus,Y, b);
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
        
        Ysub = Y1(tinx1:tinx2,:);
        
        T1 = size(Ysub,1);
        C1 = C(1:T1,1:T1);        
        [npart, status, ~]=getn(C1,a(:,1),Ysub,lamn,rel_tol,1);
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
        Yn = bsxfun(@minus,mean(Y,2),b);
        history.lamns = NaN;
        history.lambda1s = NaN;
        history.lambda2s = NaN;
        history.posts = NaN; 
        fprintf('Fail\n');
        return;
    else
        if N>1
            % L2 norm
%             P = Cn;
%             Q = P'*P+1/lambda1;
%             K = P'*Y1 ;            
%             aT = Q\K;            
%             a = aT';
%             b = mean(mean(Y-Cn*aT));
%             
            %L1-norm            
            [a,status,~] = geta(Cn,Y,1/lambda1,rel_tol,1);            
%             b=0;
        else
            a=1;
        end
        b=(mean(Y-Cn*a'));
        E = bsxfun(@minus,Y-Cn*a(:,1)',b);    
        e2 = E(:)'*E(:);    
        lamns(iter)=2*lambda2/lambda1;    
        lambda1s(iter) = lambda1;
        lambda2s(iter) = lambda2;
        as(:,iter) = a(:,1);
        bs(:,iter) =b;
        ns(:,iter)=n;
    %     post1s(iter) = T*N/2*log(lambda1) - lambda1/2*e2;
%         posts(iter) = T*N/2*log(lambda1) - lambda1/2*e2 + T*N*log(lambda2)- N*lambda2*sum(n) - 2*N*a(:)'*a(:);
        posts(iter) = T*N/2*log(lambda1) - lambda1/2*e2 + T*log(lambda2)- lambda2*sum(n) -sum(a);
    end
    
    if (iter>1 && diff(posts(iter-1:iter))<1e-3) || (iter>2 && diff(posts(iter-1:iter))<1e-3*diff(posts(1:2))) ||iter == MAX_ITER
        % renormalization of n and a
        a1 = a/sum(a);
        Yb = bsxfun(@minus,Y,b); 
        scale = Cn'*Yb*a1/(Cn'*Cn);
        n = n*scale;
        Yn = Yb*a1;
        
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
  
    
    
    
    fprintf('.');
end




