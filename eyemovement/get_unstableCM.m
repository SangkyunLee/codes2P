function out = get_unstableCM(CM,ich,THR)
% THR(1): eyemove/frame
% THR(2): eyemove/frame, 
% default THR=[20 5]
if nargin<3
    THR=[20 5 40];
end

df_PARs=[0 diff(CM(ich,:));diff(CM(ich,:))  0];
dist = sqrt(sum(CM(1:2,:).^2,1));
rdist =dist -nanmean(dist);


inxP1 = (abs(df_PARs(1,:)>THR(2) & (df_PARs(1,:).*df_PARs(2,:))<0));
inxP2 = abs(df_PARs(1,:))>THR(1) ;
inxP = inxP1 | inxP2;
inxbad = find( inxP | rdist>THR(3) );
inxbad = setdiff(inxbad, [1 size(CM,2)]);


nh = [inxbad-1; inxbad+1];
out = inxbad;
for i= 1 : size(nh,1)
    D = CM(ich,inxbad)-CM(ich,nh(i,:));
    out =[out nh(i,abs(D)>10)];
end
out = unique(out);