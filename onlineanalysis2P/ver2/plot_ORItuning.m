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
Nrep = Ntrials /Ncond;
stimtime = timeinfo.stimtime;
frame_start = timeinfo.frame_start;


 
 
 %% plot average response
 imgfcond = stimtime(frame_start>0); % stimulus condition of image frames
uevt = unique(others.events{1});
nevt = length(uevt);
Y = permute(Y{1},[2 1 3]);
nf = size(Y,1); % number of frame
nC = size(Y,3); % number of cell
AY = zeros(nf,nevt,nC);
for i = 1: nevt
    inx = others.events{1} == uevt(i);
    AY(:,i,:)=squeeze(mean(Y(:,inx,:),2));
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
plot_condResptime(AY,ACond, tpar)


%% plot mean response across orientations
AY2 = squeeze(mean(Y,2));
tpar2 = tpar;
tpar2.imgfcond = imgfcond2(1:size(AY2,1));
tpar2.imgfrep = ones(size(AY2,1),1);
ACond2.Ncond=1;
ACond2.Ntrials =1;
plot_condResptime(AY2,ACond2, tpar2)

%% plot tuning function
NX=2; % xgrid 
NY=3; % ygrid
% ntr = 5; % 
sframe = 5:9; % frames selected for ori tuning function

Resp = zeros(nevt,nC);
EResp = zeros(nevt,nC);

for i = 1: nevt
    inx = others.events{1} == uevt(i);
    Resp_tr = mean(Y(sframe,inx,:),1);
    Resp(i,:)= squeeze(mean(Resp_tr,2));
    EResp(i,:) = squeeze(std(Resp_tr,0,2))/sqrt(Nrep);
end
ORI = Params.VStimparam.orientation;
figure; 
for ic = 1: nC
    subplot(NY,NX,ic);
    plot(ORI,Resp(:,ic));
end

% 
% 
%  
