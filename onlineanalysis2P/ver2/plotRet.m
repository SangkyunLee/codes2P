data = struct;
data.dFF = dFF;
data.timeinfo = timeinfo;
x=round([-0.5 4]/Params.msperframe*1000);
spec.frames = x(1):x(end);
spec.dataType ='dFF';
spec.nCell=  size(dFF,2);
[Y, others]=data_sort(data,spec,0, Params);




%-------- get stimulus condition and repetition based on image frames
Ntrials = length(Cond.condseq);
Ncond = Cond.Ncond;
% Nx = Cond.Nx;
Nrep = Ntrials /Ncond;
% Ny = Ncond/Nx;
stimtime = timeinfo.stimtime;
frame_start = timeinfo.frame_start;


%% plot signal across time
imgfcond = stimtime(frame_start>0); % stimulus condition of image frames
rep = get_stimrep_imgf(stimtime,Nrep,Ncond);
imgfrep =  rep(frame_start>0);
tpar.imgfcond = imgfcond;
tpar.imgfrep = imgfrep;
tpar.msperframe = Params.msperframe;
 plot_onRettime(dFF,Cond, tpar)
 
 
 %% plot average response
 
uevt = unique(others.events{1});
nevt = length(uevt);
Y = Y{1};
nf = size(Y,2);
nC = size(Y,3);
AY = zeros(nf,nevt,nC);
for i = 1: nevt
    inx = others.events{1} == uevt(i);
    AY(:,i,:)=squeeze(mean(Y(inx,:,:),1));
end
AY = reshape(AY, [nf*nevt, nC]);
ACond = Cond;
ACond.Ntrials = Ncond;


durf = round(length(imgfcond(imgfcond==uevt(1)))/Nrep);
inxst = find(spec.frames==0);
inx2 = inxst : inxst+durf-1;
imgfcond2 = zeros(nf,nevt);
imgfcond2(inx2,:) = repmat(uevt(:)',[length(inx2) 1]);
imgfcond2 = imgfcond2(:);

tpar.imgfcond = imgfcond2;
tpar.imgfrep = ones(size(imgfcond2));
tpar.msperframe = Params.msperframe;
plot_onRettime(AY,ACond, tpar)


 
