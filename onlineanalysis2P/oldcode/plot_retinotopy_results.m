function plot_retinotopy_results(dFF,Params,Cond,plotopt)
% function plot_retinotopy_results(dFF,Params,Cond,plotopt)
% fix a bug (incorrect time info), 2014-01-26, Sangkyun Lee

Ntrials = Cond.Ntrials;
Ncond = Cond.Ncond;
Nx = Cond.Nx;
% Ny = Cond.Ny;


%% estimation of the spike rates
if isfield(plotopt,'AP') 
    if plotopt.AP==1,
    
        Nhats=[]; Phats=[]; Chats=[];
        for iroi=1:size(dFF,2)
            F = dFF(:,iroi);
            clear V C
            T       = size(F,1); % # of time steps
            V.dt    = Params.msperframe/1000;  % time step size
            V.Ncells = size(F,2);
            V.Npixels = size(F,2);
            [Nhat Phat Vhat Chat]=fast_oopsi(F,V);
            Nhats=[Nhats Nhat];
            Phats=[Phats Phat];
            Chats=[Chats Chat];
        end


%         normF=(F)/(max(F));
%         normC= (Chat)/(max(Chat));
%         sig=[normF normC Nhat/max(Nhat)];
%         offsets_ch = [size(sig,2)-1:-1:0]/2;
%         sig = sig + repmat(offsets_ch,size(sig,1),1);
        
        sig = Nhats;
    end
else
    sig =dFF;
end
        

%%
a = zeros(size(Params.timeNI));
for ii = 1: length(Params.stim_onset_tstampinNI)    
    inxst = Params.stim_onset_tstampinNI(ii);
    inxed = inxst + Params.stim_duration_tstampinNI(ii)-1;
    a(inxst:inxed)=ii;  
end

b = zeros(Params.Nframes,2);
imagframe_stamp = Params.mirror_start_time;
imagdur = Params.diff_mirror_start_time;
if length(imagdur)==Params.Nframes-1,
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


Nrepetition = Ntrials /Ncond;
maxV= max(sig(:));
minV = min(sig(:));
for ii=1:Nx
    clr=0.4+(Nx-ii)*ones(1,3)*0.5/Nx;
    colors{ii} = clr;
    colors{ii+Nx} = clr;
    colors{ii+2*Nx} = clr;
end
colors{Ncond} = [0.1 0.1 0.1];
frametimes = Params.frame_times;

figure('Color','white')
for ic =1 : Ntrials
    inxs=find(condtimeseries==ic);    
    X=[inxs', inxs(end:-1:1)']*Params.msperframe/1000;
    Y=[maxV*ones(1,length(inxs)) minV*ones(1,length(inxs))];
    icolor = mod(ic,Ncond);
    if icolor ==0, icolor = Ncond; end
    set(fill(X,Y,colors{icolor}),'EdgeColor',colors{icolor});
    hold on;    
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
figure('Color','white')
accumN =[];
for ic =1 : Ncond
        
    for iter=1:Nrepetition
        inxs=find(condtimeseries==ic+Ncond*(iter-1));
        subsig = sig(inxs,:);
        Nsubsig = size(subsig,1)
        
        if iter ==1,
            sig_conds{ic} = subsig;
            minNsubsig = Nsubsig;
        else
            if Nsubsig < minNsubsig,            
                sig_conds{ic} = sig_conds{ic}(1:Nsubsig,:) + subsig;
                minNsubsig = Nsubsig;
            else
                sig_conds{ic} = sig_conds{ic} + subsig(1:minNsubsig,:);
            end
        end        
    end   

    X=[[1:minNsubsig], [minNsubsig:-1:1]]+sum(accumN);
    X = X*Params.msperframe/1000;
    Y=[maxV*ones(1,minNsubsig) minV*ones(1,minNsubsig)];
    accumN = [ accumN minNsubsig];
    icolor = mod(ic,Ncond);
    if icolor ==0, icolor=Ncond; end
    set(fill(X,Y,colors{icolor}),'EdgeColor',colors{icolor});
    hold on;  
    plot(X(1:end/2),sig_conds{ic}/Nrepetition);
    
end

set(gca, 'YTick',fliplr(offsets_ch));
Yticklabel={};
for ich=1:size(sig,2)
    Yticklabel{ich} = num2str(size(sig,2)-ich+1);
end
 set(gca,'YTickLabel',Yticklabel)
 minV=min(sig_conds{1}(:)/Nrepetition);
maxV=max(sig_conds{1}(:)/Nrepetition);
ylim([0 1.1*maxV]);
% xlim([0 max(X)])
xlabel('Time(sec)','FontSize',14)
ylabel('Cell','FontSize',14)
set(gca,'FontSize',14);