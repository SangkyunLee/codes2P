function h = tent(T, inth)
% function h = tent(T, inth)
% generate picewise linear basis functions
% T: number of time samples
% inth: interval of the half tent
% Sangkyun Lee


xgrd = [1:T];
xcs = [1:inth:T];
h = zeros(T,length(xcs));

kk=1;
% figure;
for xc=xcs
    
    hi = zeros(1,T);
    x1 = (xc)-inth;
    fxi = (xgrd-x1)/inth;
    mask1 =[x1:x1+inth];
    inxm=find( mask1>0 & mask1<=T)  ;  
    hi(mask1(inxm))=fxi(mask1(inxm));
    
    x2 = xc+inth;
    fxi = (x2-xgrd)/inth;
    mask2 =[xc:x2];
    inxm=find( mask2>0 & mask2<=T) ;   
    hi(mask2(inxm))=fxi(mask2(inxm));
    h(:,kk)=hi(:);
    kk = kk+1;
%      plot(hi); hold on;
end

