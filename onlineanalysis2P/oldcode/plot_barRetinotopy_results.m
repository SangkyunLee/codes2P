function plot_barRretinotopy_results(dFF,Params)

sig =dFF;
stimparam = Params.stimparam;
nstepin1dir = size(Params.stim_onset_tstampinNI,1);
nrep = stimparam.nrepetition;
ndir = length(stimparam.dir);  
if isfield(stimparam,'blank_samplesinNI') && stimparam.blank_samplesinNI>0    
    stim_onset_tstampinNI =  Params.stim_onset_tstampinNI(1:end-1,:);
    nstepin1dir = nstepin1dir-1;
end
%%
a = zeros(size(Params.timeNI));
for ii = 1: nrep
    for jj = 1 :ndir
        for kk = 1: nstepin1dir    
            inxst = stim_onset_tstampinNI(kk,(ii-1)*ndir+jj);
            inxed = inxst + Params.stim_duration_tstampinNI(kk,(ii-1)*ndir+jj)-1;
            a(inxst:inxed) = (jj-1)*nstepin1dir + kk; 
        end
    end
end

b = zeros(Params.Nframes,2);
imagframe_stamp = Params.mirror_start_time;
imagdur = Params.diff_mirror_start_time;
if imagdur==Params.Nframes-1,
    imagdur(end+1)=mean(imagdur);
end
inxmulticonds = [];
for ii = 1: Params.Nframes
    inxst = imagframe_stamp(ii);
    inxed = inxst + imagdur(ii)-1;
    listconds = a(inxst:inxed);
    conds=unique(listconds);
    if length(conds)>1, 
        inxdiffcond = find(diff(listconds)~=0);
        ratio = round(inxdiffcond/(inxed-inxst+1)*100); 
        inxmulticonds=[inxmulticonds; [ii listconds(inxdiffcond) listconds(inxdiffcond+1) ratio ]]; 
    end
    b(ii,1:length(conds))=conds;    
end


condtimeseries =b(:,1);
for im=1: size(inxmulticonds,1)
    if inxmulticonds(im,4)>50,
        condtimeseries(inxmulticonds(im,1)) =  inxmulticonds(im,2);
    else
        condtimeseries(inxmulticonds(im,1)) =  inxmulticonds(im,3);
    end
end



%% plot signal across time
% sig = dFF;
offsets_ch = [size(sig,2)-1:-1:0]/2;
sig = sig + repmat(offsets_ch,size(sig,1),1);

Ncond = nstepin1dir*ndir
maxV= max(sig(:));
minV = min(sig(:));
for ii=1:nstepin1dir
    clr=0.4+(nstepin1dir-ii)*ones(1,3)*0.5/nstepin1dir;
    for jj = 1 : ndir
        colors{ii+(jj-1)*nstepin1dir} = clr;
    end
end
colors{Ncond} = [0.1 0.1 0.1];
frametimes = Params.frame_times;

figure('Color','white')
for ic =1 : Ncond 
    inxs=find(condtimeseries==ic);    

    stinx = [inxs(1); inxs(find(diff(inxs)>1)+1)];
    edinx = [inxs(find(diff(inxs)>1)); inxs(end)];
    
    for ic2 = 1 : length(stinx)
        inx1 = stinx(ic2):edinx(ic2);
        X=[inx1 inx1(end:-1:1)]*Params.msperframe/1000;
        Y=[maxV*ones(1,length(inx1)) minV*ones(1,length(inx1))];
        icolor = mod(ic,Ncond)
        if icolor ==0, icolor = Ncond; end
        set(fill(X,Y,colors{icolor}),'EdgeColor',colors{icolor});
        hold on;    
    end
end
plot((frametimes(1:Params.Nframes)-frametimes(1))/Params.samplingfreq_NI,sig);

set(gca, 'YTick',fliplr(offsets_ch));
Yticklabel={};
for ich=1:size(sig,2)
    Yticklabel{ich} = num2str(size(sig,2)-ich+1);
end
 set(gca,'YTickLabel',Yticklabel)
ylim([minV maxV]);xlim([0 frametimes(Params.Nframes)/Params.samplingfreq_NI])
xlabel('Time(sec)','FontSize',14)
ylabel('Cell','FontSize',14)
set(gca,'FontSize',14);


%% plot dF/F across conditions

df = diff(condtimeseries);
inx1 = find(condtimeseries==1);
inx2  = find(df==min(df));

for ii=1:length(inx2)
    if ii==1,
        inx_iter(ii,1) = inx1(1);
    else        
        inx3 = min(abs(inx1 - inx2(ii-1)));        
        inx_iter(ii,1) = inx2(ii-1) + inx3; 
    end
    if ii==length(inx2)
        mblankinx = round(mean(blankinxs));
        inx_iter(ii,2) = inx2(ii)+mblankinx -1;
    else
        blankinx = min(abs(inx1 - inx2(ii)));
        blankinxs(ii)=blankinx;
        inx_iter(ii,2) = inx2(ii)+blankinx -1
    end
end



    

        
for iter=1:nrep
    
    subsig = sig(inx_iter(iter,1):inx_iter(iter,2),:);
    Nsubsig = size(subsig,1)

    if iter ==1,
        sigiter = subsig;
        minNsubsig = Nsubsig;
    else
        if Nsubsig < minNsubsig,            
            sigiter = sigiter(1:Nsubsig,:) + subsig;
            minNsubsig = Nsubsig;
        else
            sigiter = sigiter + subsig(1:minNsubsig,:);
        end
    end        
end   

figure('Color','white')
for ic =1 : Ncond 
    inxs=find(condtimeseries(inx_iter(1,1):minNsubsig + inx_iter(1,1)-1)==ic);    

    stinx = [inxs(1); inxs(find(diff(inxs)>1)+1)];
    edinx = [inxs(find(diff(inxs)>1)); inxs(end)];
    
    for ic2 = 1 : length(stinx)
        inx1 = stinx(ic2):edinx(ic2);
        X=[inx1 inx1(end:-1:1)]*Params.msperframe/1000;
        Y=[maxV*ones(1,length(inx1)) minV*ones(1,length(inx1))];
        icolor = mod(ic,Ncond)
        if icolor ==0, icolor = Ncond; end
        set(fill(X,Y,colors{icolor}),'EdgeColor',colors{icolor});
        hold on;    
    end
end

plot((1:minNsubsig)*Params.msperframe/1000,sigiter/5);
    


set(gca, 'YTick',fliplr(offsets_ch));
Yticklabel={};
for ich=1:size(sig,2)
    Yticklabel{ich} = num2str(size(sig,2)-ich+1);
end
 set(gca,'YTickLabel',Yticklabel)
 minV=min(sigiter(:)/nrep);
maxV=max(sigiter(:)/nrep);
ylim([0 1.1*maxV]);
% xlim([0 max(X)])
xlabel('Time(sec)','FontSize',14)
ylabel('Cell','FontSize',14)
set(gca,'FontSize',14);