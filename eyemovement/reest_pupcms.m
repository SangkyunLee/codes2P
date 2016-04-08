function [out, inxfails,outtmp] = reest_pupcms(out,data,CM0,inxfail,thrs,Niter,PL)

chunksize=1000;
% thrs =[3 5 10 15 20 30];
% Niter =1;
htwin=10;
inxfails = cell(Niter,1);
Nframe = size(data,3);

for iter=1:Niter
    CMs=zeros(2,Nframe, length(thrs));
    nzCM=zeros(Nframe, length(thrs));
    opt.winsize=10;
    opt.bdisp=false;
    for ithr = length(thrs):-1:1    
        opt.thr = thrs(ithr);
        %track_rawdata2(MOV,opts, newCM,htwin, frames)
        if PL
            hf = @(data,CM,selfr)track_rawdata2_PL(data,opt,CM,6,selfr);
        else
            hf = @(data,CM,selfr)track_rawdata2(data,opt,CM,6,selfr);
        end
        
        
        outtmp(ithr) = apply_fun_chunk (hf, chunksize,out,data,inxfail, CM0);

        CMs(:,:,ithr)=outtmp(ithr).CM;
        nzCM(:,ithr) = outtmp(ithr).CM(1,:)>0;
    end


    [newCM1, inxout] = slow_xy(CMs,htwin,CM0,inxfail,1);
    out.CM(:,inxfail)=newCM1;

    inxoutm = min(inxout,[],1);
    for ithr =1 : length(thrs)    
        opt.thr = thrs(ithr);
        ix1 = find(inxoutm==ithr);    
        out.sPIX(inxfail(ix1)) = outtmp(ithr).sPIX(inxfail(ix1));
    end  
    
    
    inxfail = get_unstableCM(out.CM,1);
    inxfails{iter} = inxfail;
end