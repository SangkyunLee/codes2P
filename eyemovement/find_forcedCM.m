function out = find_forcedCM(CM,selfr,CM0)

M= squeeze(CM(1,selfr,:));
M0 = M-CM0(1,1);
[~,r] = min(abs(M0),[],2);
inx =sub2ind(size(M),1:size(M,1),r');
CMx = M(inx);
CMx(CMx==0) =CM0(1,1);


M= squeeze(CM(2,selfr,:));
M0 = M-CM0(2,1);
[~,r] = min(abs(M0),[],2);
inx =sub2ind(size(M),1:size(M,1),r');
CMy = M(inx);
CMy(CMy==0) =CM0(2,1);

out =[CMx; CMy];