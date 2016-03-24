function [X0 Y0]=cmass(L)

X_hist=sum(L,1); 
Y_hist=sum(L,2); 
X=1:size(L,2); Y=1:size(L,1); 
X0 = sum(X.*X_hist)/sum(X_hist); 
Y0 = sum(Y'.*Y_hist)/sum(Y_hist);