function [W out]= sl_smlr_matlab(W,X,Y,Bkkc,delta_kkc,lambda1, lambda2,MAX_ITER)

[m,n] = size(Y);[d,n] = size(X);

 WtX=W'*X;  %[d x m-1]
    denom = sum(exp(WtX));
    LL=-Inf*ones(1,MAX_ITER);
    Lp=-Inf*ones(1,MAX_ITER);
    Wdiff=-Inf*ones(1,MAX_ITER);

    for iter =1: MAX_ITER      
    % if mod(iter,100)==0, iter
    % end

        Wdiff(iter)=0;
        for inx_i = 1 : m-1 % number of weight vectors        
            for inx_j = 1 : d % input vector dimension
                if inx_j == 1,
                    WiX=W(:,inx_i)'*X;
                    exp_WiX = exp(WiX); 
                    Prob_i = exp_WiX./denom;
                end


                Ai = Y(inx_i,:)-Prob_i;
                grad_Wi = X*Ai';

                Wji_new = (Bkkc(inx_j)*W(inx_j,inx_i) - grad_Wi(inx_j))/(Bkkc(inx_j)-lambda1);
                
                Wji_new = softmax(Wji_new,delta_kkc(inx_j));
                
                
                if isnan(Wji_new)
                    keyboard;
                end
                Wji_old = W(inx_j,inx_i);            
                
                Wdiff(iter) = Wdiff(iter)+ abs(Wji_new-Wji_old)/((m-1)*d);

    %             if abs(Wji_old -Wji_new)<1e-10,
    %                 a=[inx_i, inx_j abs(Wji_old -Wji_new)]
    %             end


                WiX_new = WiX + (Wji_new - Wji_old)*X(inx_j,:);  %%WiX_new = W(:,inx_i)'*X;
                exp_WiX_new = exp(WiX_new);
                denom = denom - exp_WiX + exp_WiX_new;  %%denom = sum(exp(WtX)); <--WtX=W'*X;

                % update
                W(inx_j,inx_i)=Wji_new;
                Prob_i = exp_WiX_new./denom;            
                WiX = WiX_new;
                exp_WiX = exp_WiX_new; 

    %             W(end-5:end,:)
            end
        end 
        b=W'*X;
        aa=sum(exp(b));
        LL(iter)= sum(sum(Y.*(b))) - sum(log(aa));
        Lp(iter) = LL(iter) - lambda1*sum(sum(W.^2)) - lambda2*sum(sum(abs(W)));
        
        if iter>1,
            diff=abs(Lp(iter)-Lp(iter-1));
            if iter==2, diff1=diff; end
            if diff<1e-3 || (diff/diff1)<1e-5,
                disp(['Converged in iteration:' num2str(iter)]);                
                break;
            end
        end                
        
    end
    
    out=[LL; Lp; Wdiff];