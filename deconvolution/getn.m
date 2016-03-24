function [x,status,history] = getn(C,a,Y,lambda,varargin)
%
% function [x,status,history] = getn(C,a,Y,lambda,varargin)
%
%       minimize   ||C*x*a'-Y||^2 + lambda*sum(x_i),
%       subject to x_i >= 0, i=1,...,n 
%
% INPUT
%   C       : TxT matrix; Calcium convolution matrix.
%   a       : m vector; scale of each pixel
%
%   Y       : T x m Fluorescence matrix.
%   lambda  : positive scalar; regularization parameter
%
%   tar_gap : relative target duality gap (default: 1e-3)
%   quiet   : boolean; suppress printing message when true (default: false)
% OUTPUT
%   x       : T vector; spike information
%   status  : string; 'Solved' or 'Failed'
%
%   history : matrix of history data. columns represent (truncated) Newton
%             iterations; rows represent the following:
%            - 1st row) gap
%            - 2nd row) primal objective
%            - 3rd row) dual objective
%            - 4th row) step size
%            - 5th row) pcg iterations
%            - 6th row) pcg status flag
%
%
% This code is a modifed version of
% l1-Regularized Least Squares Problem Solver.
% The original version of the l1-least square code was written
% by Kwangmoo Koh, Seung-Jean Kim, and Stephen Boyd
%
% 10-23-2013, modified by Sangkyun Lee
%
% IPM PARAMETERS
MU              = 2;        % updating parameter of t
MAX_NT_ITER     = 400;      % maximum IPM (Newton) iteration

% LINE SEARCH PARAMETERS
ALPHA           = 0.01;     % minimum fraction of decrease in the objective
BETA            = 0.5;      % stepsize decrease factor
MAX_LS_ITER     = 100;      % maximum backtracking line search iteration


    

Ct = C';
[~,T] = size(C);
% if T~=n, 
%     error('matrix C should be square');
% end
if size(a,1)~=size(Y,2)
    error('The size of vector a should be the same as the number of columns of Y');
end
% VARIABLE ARGUMENT HANDLING
t0         = min(max(1,1/lambda),T/1e-3);
defaults   = {1e-3,false,1e-3,5000,ones(T,1),t0};
given_args = ~cellfun('isempty',varargin);
defaults(given_args) = varargin(given_args);
[reltol,quiet,eta,pcgmaxi,x,t] = deal(defaults{:});

f = -x;

% RESULT/HISTORY VARIABLES
pobjs = zeros(1,MAX_NT_ITER+1); 
dobjs = zeros(1,MAX_NT_ITER+1); 
sts = zeros(1,MAX_NT_ITER+1); 
pitrs = zeros(1,MAX_NT_ITER+1); 
pflgs = zeros(1,MAX_NT_ITER+1);


pobj  = Inf; dobj  =-Inf; s   = Inf; pitr  = 0 ; pflg  = 0 ;

ntiter  = 0; lsiter  = 0; 
dx =  zeros(T,1);

% diagxtx = diag(Ct*C);
diagxtx = 2*ones(T,1);

if (~quiet) 
    fprintf('\nSolving a problem of size (T=%d), with lambda=%.5e\n',T,lambda);
    fprintf('-----------------------------------------------------------------------------\n');
    fprintf('%5s %9s %15s %15s %13s %11s\n',...
            'iter','gap','primobj','dualobj','step len','pcg iters'); 
end

%------------------------------------------------------------
%               MAIN LOOP
%------------------------------------------------------------

for ntiter = 0:MAX_NT_ITER
    
    Z = C*x*a'-Y;
    
    %------------------------------------------------------------
    %       CALCULATE DUALITY GAP
    %------------------------------------------------------------


    Nu = 2*Z;
    CNua = C'*Nu*a;    
    minCnu = min(CNua);    
    
    if (minCnu < -lambda)
        Nu = Nu*lambda/(-minCnu);
    end
    
    pobj  =  Z(:)'*Z(:) + lambda*sum(x(1:end-1),1);
    dobj  =  max(-0.25*Nu(:)'*Nu(:)-Nu(:)'*Y(:),dobj);
    gap   =  pobj - dobj;

    pobjs(ntiter+1) = pobj; dobjs(ntiter+1) =  dobj; sts(ntiter+1) = s;
    pflgs(ntiter+1) = pflg; pitrs(ntiter+1) = pitr;


    %------------------------------------------------------------
    %   STOPPING CRITERION
    %------------------------------------------------------------
    if (~quiet)
        fprintf('%4d %12.2e %15.5e %15.5e %11.1e %8d\n',...
        ntiter, gap, pobj, dobj, s, pitr);
    end

    if (gap/abs(dobj) < reltol) 
        status  = 'Solved';
        pobjs = pobjs(1:ntiter+1); dobjs = dobjs(1:ntiter+1); sts = sts(1:ntiter+1);
        pflgs = pflgs(1:ntiter+1); pitrs = pitrs(1:ntiter+1);
        history = [pobjs-dobjs; pobjs; dobjs; sts; pitrs; pflgs];
        if (~quiet) 
            disp('Absolute tolerance reached.');
        end
        %disp(sprintf('total pcg iters = %d\n',sum(pitrs)));
        return;
    end
    %------------------------------------------------------------
    %       UPDATE t
    %------------------------------------------------------------
    if (s >= 0.5)
        t = max(min(T*MU/gap, MU*t), t);
    end

    %------------------------------------------------------------
    %       CALCULATE NEWTON STEP
    %------------------------------------------------------------
    
    d1 = (1/t)./(x.^2);
%     d1 = [(1/t)./(x(1:end-1).^2); 0];

    % calculate gradient
    gradphi = Ct*(Z*2)*a + lambda-(1/t)./x;
%     gradphi = Ct*(Z*2)*a + [lambda-(1/t)./x(1:end-1); 0];
    
    % calculate vectors to be used in the preconditioner
    prb     = diagxtx+d1;

    % set pcg tolerance (relative)
    normg   = norm(gradphi);
    pcgtol  = min(1e-1,eta*gap/min(1,normg));
    
    if (ntiter ~= 0 && pitr == 0) 
        pcgtol = pcgtol*0.1; 
    end


    [dx,pflg,~,pitr,~] = ...
        pcg(@AXfunc_l1_ls,-gradphi,pcgtol,pcgmaxi,@Mfunc_l1_ls,...
            [],dx,C,Ct,a,d1,1./prb);

    %dx = (2*C'*C+diag(d1))\(-gradphi);

    if (pflg == 1) 
        pitr = pcgmaxi; 
    end
    
    %------------------------------------------------------------
    %   BACKTRACKING LINE SEARCH
    %------------------------------------------------------------
    phi = Z(:)'*Z(:)+lambda*sum(x)-sum(log(-f))/t;
    s = 1.0;
    gdx = gradphi'*dx;
    for lsiter = 1:MAX_LS_ITER
        newx = x+s*dx;
        newf = -newx;
        if (max(newf) < 0)
            newZ   =  C*newx*a'-Y;
            newphi =  newZ(:)'*newZ(:) + lambda*sum(newx)-sum(log(-newf))/t;
            if (newphi-phi <= ALPHA*s*gdx)
                break;
            end
        end
        s = BETA*s;
    end
    if (lsiter == MAX_LS_ITER) 
        break;
    end % exit by BLS
        
    x = newx; f = newf;
end


%------------------------------------------------------------
%       ABNORMAL TERMINATION (FALL THROUGH)
%------------------------------------------------------------
if (lsiter == MAX_LS_ITER)
    % failed in backtracking linesearch.
    if (~quiet) 
        disp('MAX_LS_ITER exceeded in BLS');
    end
    status = 'Failed';
elseif (ntiter == MAX_NT_ITER)
    % fail to find the solution within MAX_NT_ITER
    if (~quiet)
        disp('MAX_NT_ITER exceeded.'); 
    end
    status = 'Failed';
end
history = [pobjs-dobjs; pobjs; dobjs; sts; pitrs; pflgs];

return;

%------------------------------------------------------------
%       COMPUTE AX (PCG)
%------------------------------------------------------------
function [Y] = AXfunc_l1_ls(x,C,Ct,a,d1,p1)

Y = (Ct*((C*x)*2)*a'*a)+d1.*x;

%------------------------------------------------------------
%       COMPUTE P^{-1}X (PCG)
%------------------------------------------------------------
function [Y] = Mfunc_l1_ls(x,C,Ct,a,d1,p1)

Y = p1.*x;
