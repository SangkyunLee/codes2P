function plot_onRettime(sig,Cond, tpar)
    
    imgfcond = tpar.imgfcond;
    imgfrep = tpar.imgfrep;
    msperframe = tpar.msperframe;
    Nframes = length(imgfcond);

    Ncond = Cond.Ncond;
    Ntrials = length(Cond.condseq);
    Nx = Cond.Nx;
    Ny = Ncond/Nx;
    Nrep = Ntrials/Ncond;


    
    
    offsets_ch = (size(sig,2)-1:-1:0)/2;
    sig = sig + repmat(offsets_ch,size(sig,1),1);

    maxV= max(sig(:));
    minV = min(sig(:));
    colors = cell(Ncond,1);
    for ii=1:Nx
        clr=0.4+(Nx-ii)*ones(1,3)*0.5/Nx;
        for j= 0 : Ny-1
            colors{ii +j*Nx} = clr;
        end
    end
    colors{Ncond} = [0.1 0.1 0.1];


    figure('Color','white')
    for irep = 1: Nrep
        for ic =1 : Ncond
            inxs=find(imgfcond==ic & imgfrep ==irep);    
            X=[inxs', inxs(end:-1:1)']*msperframe/1000;
            Y=[maxV*ones(1,length(inxs)) minV*ones(1,length(inxs))];
            icolor = mod(ic,Ncond);
            if icolor ==0, icolor = Ncond; end
            set(fill(X,Y,colors{icolor}),'EdgeColor',colors{icolor});
            hold on;    
        end
    end
    tframe = (0:Nframes-1)*msperframe/1000;
    plot(tframe,sig);

    set(gca, 'YTick',fliplr(offsets_ch));
    Yticklabel=cell(size(sig,2),1);
    for ich=1:size(sig,2)
        Yticklabel{ich} = num2str(size(sig,2)-ich+1);
    end
     set(gca,'YTickLabel',Yticklabel)
    ylim([minV maxV]);xlim([0 tframe(Nframes)])
    xlabel('Time(sec)','FontSize',14)
    ylabel('Cell','FontSize',14)
    set(gca,'FontSize',14);

end



