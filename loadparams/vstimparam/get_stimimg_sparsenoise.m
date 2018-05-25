function X = get_stimseq_sparsenoise(vstimparam,tinfo)
% function X = get_stimseq_sparsenoise(vstimparam)
    bgcolor = 127.5;

    stimimg = get_stim2Dimg(vstimparam);
    
    Nframe = data(idata).Params.Nframes;

    frame_evt = tinfo.stimtime(tinfo.frame_start>0);

    X = bgcolor*ones(size(stimimg,1),size(stimimg,2),Nframe);

    for i = 1: Nframe
        ievt = frame_evt(i);
        if ievt>0,
            X(:,:,i) = stimimg(:,:,ievt);
        end
    end

end


function stimimg = get_stim2Dimg(vstimparam, bgcolor)
    % stimseq: Y x X x condition

    evtcond = vstimparam.evtcond;
    dotNumX = vstimparam.dotNumX;
    dotNumY = vstimparam.dotNumY;
    [a, b]=unique(evtcond(:,end));
    dotClr = cell2mat(vstimparam.dotColors);
    dotClr = unique(dotClr);

    Ncond = length(a);
    evtmap = evtcond(b,:);
    stimimg = bgcolor * ones(dotNumY,dotNumX,Ncond);

    for i= 1 : Ncond
        ix = evtmap(i,1);
        iy = evtmap(i,2);
        ic = evtmap(i,3);
        stimimg(iy,ix,i) = dotClr(ic);
    end
end

