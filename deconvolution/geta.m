function [x,status,history] = geta(Cn,Y,lambda,varargin)
%
% IPM PARAMETERS
MU              = 2;        % updating parameter of t
MAX_NT_ITER     = 400;      % maximum IPM (Newton) iteration

% LINE SEARCH PARAMETERS
ALPHA           = 0.01;     % minimum fraction of decrease in the objective
BETA            = 0.5;      % stepsize decrease factor
MAX_LS_ITER     = 100;      % maximum backtracking line search iteration


Cnt = Cn';
m = size(Cn,1);        
n = size(Y,2);
    


% VARIABLE ARGUMENT HANDLING
t0         = min(max(1,1/lambda),n/1e-3);
defaults   = {1e-3,false,1e-3,5000,ones(n,1),t0};
given_args = ~cellfun('isempty',varargin);
defaults(given_args) = varargin(given_args);
[reltol,quiet,eta,pcgmaxi,x,t] = deal(defaults{:});

f = -x;

% RESULT/HISTORY VARIABLES
pobjs = [] ; dobjs = [] ; sts = [] ; pitrs = []; pflgs = [];
pobj  = Inf; dobj  =-Inf; s   = Inf; pitr  = 0 ; pflg  = 0 ;

ntiter  = 0; lsiter  = 0; zntiter = 0; zlsiter = 0;
normg   = 0; prelres = 0; dx =  zeros(n,1);

% diagxtx = diag(Ct*C);
diagxtx = 2*ones(n,1);

if (~quiet) disp(sprintf('\nSolving a problem of size (m=%d, n=%d), with lambda=%.5e',...
            m,n,lambda)); end
if (~quiet) disp('-----------------------------------------------------------------------------');end
if (~quiet) disp(sprintf('%5s %9s %15s %15s %13s %11s',...
            'iter','gap','primobj','dualobj','step len','pcg iters')); end

%------------------------------------------------------------
%               MAIN LOOP
%------------------------------------------------------------

for ntiter = 0:MAX_NT_ITER
    
    Z = Cn*x'-Y;
    
    %------------------------------------------------------------
    %       CALCULATE DUALITY GAP
    %------------------------------------------------------------


    Nu = 2*Z;
    CnNu = Cn'*Nu;    
    minCnu = min(CnNu);    
    
    if (minCnu < -lambda)
        Nu = Nu*lambda/(-minCnu);
    end
    
    pobj  =  Z(:)'*Z(:) + lambda*sum(x,1);
    dobj  =  max(-0.25*Nu(:)'*Nu(:)-Nu(:)'*Y(:),dobj);
    gap   =  pobj - dobj;

    pobjs = [pobjs pobj]; dobjs = [dobjs dobj]; sts = [sts s];
    pflgs = [pflgs pflg]; pitrs = [pitrs pitr];

    %------------------------------------------------------------
    %   STOPPING CRITERION
    %------------------------------------------------------------
    if (~quiet) disp(sprintf('%4d %12.2e %15.5e %15.5e %11.1e %8d',...
        ntiter, gap, pobj, dobj, s, pitr)); end

    if (gap/abs(dobj) < reltol) 
        status  = 'Solved';
        history = [pobjs-dobjs; pobjs; dobjs; sts; pitrs; pflgs];
        if (~quiet) disp('Absolute tolerance reached.'); end
        %disp(sprintf('total pcg iters = %d\n',sum(pitrs)));
        return;
    end
    %------------------------------------------------------------
    %       UPDATE t
    %------------------------------------------------------------
    if (s >= 0.5)
        t = max(min(n*MU/gap, MU*t), t);
    end

    %------------------------------------------------------------
    %       CALCULATE NEWTON STEP
    %------------------------------------------------------------
    
    d1 = (1/t)./(x.^2);

    % calculate gradient
    gradphi = [Cnt*(Z*2) + lambda-(1/t)./x']';
    
    % calculate vectors to be used in the preconditioner
    prb     = diagxtx+d1;

    % set pcg tolerance (relative)
    normg   = norm(gradphi);
    pcgtol  = min(1e-1,eta*gap/min(1,normg));
    
    if (ntiter ~= 0 && pitr == 0) pcgtol = pcgtol*0.1; end

if 1
    [dx,pflg,prelres,pitr,presvec] = ...
        pcg(@AXfunc_l1_ls,-gradphi,pcgtol,pcgmaxi,@Mfunc_l1_ls,...
            [],dx,Cn,Cnt,d1,1./prb);
end
    %dx = (2*C'*C+diag(d1))\(-gradphi);

    if (pflg == 1) pitr = pcgmaxi; end
    
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
            newZ   =  Cn*newx'-Y;
            newphi =  newZ(:)'*newZ(:) + lambda*sum(newx)-sum(log(-newf))/t;
            if (newphi-phi <= ALPHA*s*gdx)
                break;
            end
        end
        s = BETA*s;
    end
    if (lsiter == MAX_LS_ITER) break; end % exit by BLS
        
    x = newx; f = newf;
end


%------------------------------------------------------------
%       ABNORMAL TERMINATION (FALL THROUGH)
%------------------------------------------------------------
if (lsiter == MAX_LS_ITER)
    % failed in backtracking linesearch.
    if (~quiet) disp('MAX_LS_ITER exceeded in BLS'); end
    status = 'Failed';
elseif (ntiter == MAX_NT_ITER)
    % fail to find the solution within MAX_NT_ITER
    if (~quiet) disp('MAX_NT_ITER exceeded.'); end
    status = 'Failed';
end
history = [pobjs-dobjs; pobjs; dobjs; sts; pitrs; pflgs];

return;

%------------------------------------------------------------
%       COMPUTE AX (PCG)
%------------------------------------------------------------
function [Y] = AXfunc_l1_ls(x,Cn,Cnt,d1,p1)

Y = ((Cnt*Cn)*x*2)+d1.*x;

%------------------------------------------------------------
%       COMPUTE P^{-1}X (PCG)
%------------------------------------------------------------
function [Y] = Mfunc_l1_ls(x,Cn,Cnt,d1,p1)

Y = p1.*x;
