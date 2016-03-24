%function get_ICAfilter(X0,npc)

% npc = 2%length(lam2)
Q = X0*X0';
[V, lam2] = eigs(Q, rank(Q)-1);
lam2 = diag(lam2);

whitening  = V(:,1:npc)*diag(1./sqrt(lam2(1:npc)));
dewhitening = diag(sqrt(lam2(1:npc)))*V(:,1:npc)';
Xw = whitening'*X0;

    

numSamples = size(Xw,2);
maxrounds=100       
termtol = 1e-6;
B = rand(npc);
BOld = zeros(size(B));

iternum = 0;
minAbsCos = 0;

errvec = zeros(maxrounds,1);
while (iternum < maxrounds) && ((1 - minAbsCos)>termtol)
    iternum = iternum + 1;

    % Symmetric orthogonalization.
    B = (Xw * ((Xw' * B) .^ 2)) / numSamples;
    B = B * real(inv(B' * B)^(1/2));

    % Test for termination condition.
    minAbsCos = min(abs(diag(B' * BOld)));

    BOld = B;
    errvec(iternum) = (1 - minAbsCos);
end

W = B';
icasig = W*Xw;
%icasig2 = W*whitening'*X0;
W2 = W*whitening';
A = dewhitening'*B;

    
    
%% infomax
%         Xw = whitening'*X0;
%         W = eye(npc);
%         N=size(Xw,2)
%         eta=1
%         for iter=1:100
%             % Get estimated source signals, y.
%             y = Xw'*W; % wt vec in col of W.	
%             % Get estimated maximum entropy signals Y=cdf(y).
%             Y = tanh(y);
%             % Find value of function h.	
%             % h = log(abs(det(W))) + sum( log(eps+1-Y(:).^2) )/N;
%             detW = abs(det(W));
%             h = ( (1/N)*sum(sum(Y)) + 0.5*log(detW) );
%             % Find matrix of gradients @h/@W_ji ...
%             g = inv(W') - (2/N)*Xw*Y;
%             % Update W to increase h ... 
%             W = W + eta*g;
%             % Record h and magnitude of gradient ...
%             hs(iter)=h; gs(iter)=norm(g(:));
%         end;
%         W = W';
%         B = W';
%        A = dewhitening'*B;
%        icasig = -W*Xw;
 
    