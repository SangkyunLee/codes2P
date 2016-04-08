function [inxbadfit, PARs] =get_badfit(fitInfo,Psize)
% Psize(1): smallest pupilsize
% Psize(2): largest pupilsize
% Psize(3): pupilsize change/frame
% Psize(4): pupilsize change2/frame, 
% default Psize=[10 80 20 5]
if nargin<2
    Psize=[10 80 20 5];
end
Nframe = length(fitInfo);
PARs= zeros(3,Nframe);
list=zeros(Nframe,1);
for ii=1:Nframe
    if isempty(fitInfo(ii).par)
        list(ii)=1;
    else
        PARs(:,ii)=fitInfo(ii).par;
    end
end
df_PARs=[0 diff(PARs(3,:));diff(PARs(3,:))  0];
dist = sqrt(sum(PARs(1:2,:).^2,1));
rdist =dist -nanmean(dist);


inxP1 = (abs(df_PARs(1,:)>Psize(4) & (df_PARs(1,:).*df_PARs(2,:))<0));
inxP2 = abs(df_PARs(1,:))>Psize(3) ;
inxP = inxP1 | inxP2;
inxbadfit = find((PARs(3,:)<Psize(1) | PARs(3,:)>Psize(2)) | inxP | rdist>2*std(rdist) );