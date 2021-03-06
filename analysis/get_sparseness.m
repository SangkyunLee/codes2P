function s=get_sparseness(X, dir)

if dir ==2,
    X= X';
elseif dir>2
    error('This function is defined only in dir =1 or =2');
end

X2sum = sum(X.^2,1);
Xsum2 = sum(X,1).^2;
N = size(X,1);
s=(1 - 1/N*(Xsum2./X2sum))/(1-1/N);