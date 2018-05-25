function X = get_stimseq_sparsenoise(vstimparam,tinfo,sampf, tdinsec)
% function X = get_stimseq_sparsenoise(vstimparam,tinfo,sampf, tdinsec)
    
    bgcolor = 127.5;
    if nargin<3
        tdinsec = 0;
    end
    stimimg = get_stim2Dimg(vstimparam,bgcolor);
    Nframe = length(find(tinfo.frame_start>0));
    
    tinxf = find(tinfo.frame_start>0);
    % select stimulus for imageframe prior to -tdinsec
    tdframe = tinxf + round(sampf*tdinsec); 
    tdframe(tdframe<0) = 1;     
    frame_evt = tinfo.stimtime(tdframe);
    % when the first few frames are assigned to negative time,
    % all the stimulus condition to default (0)
    if tdinsec<0
        frame_evt(tdframe == 1) = 0;
    end
    
    X = bgcolor*ones(size(stimimg,1),size(stimimg,2),Nframe);

    for i = 1: Nframe
        ievt = frame_evt(i);
        if ievt>0,
            X(:,:,i) = stimimg(:,:,ievt);
        end
    end

end


function stimimg = get_stim2Dimg(vstimparam, bgcolor)
% function stimimg = get_stim2Dimg(vstimparam, bgcolor)
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

