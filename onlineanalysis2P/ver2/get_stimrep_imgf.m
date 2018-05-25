function rep = get_stimrep_imgf(stimtime,Nrep,Ncond)
% function rep = get_stimrep_imgf(stimtime,Nrep,Ncond)
    % identify the repetition change
    inx1 = find(stimtime ==Ncond);
    inx2 = diff(diff(inx1))>10;
    inxed = inx1(inx2)+1;
    lastrepstim =stimtime(inxed(end)+1:end);
    inxed(Nrep) = inxed(end)+find(lastrepstim==max(lastrepstim),1,'last');

    rep = zeros(size(stimtime)); % stim repetition information in NI
    inxst = zeros(Nrep,1);
    for i = 1 : Nrep
        if i==1, ist =1;
        else ist = inxed(i-1)+1; end

        ied = inxed(i);
        if i==1,
            inxst(i) = find(stimtime(ist:ied)==1,1,'first');
        else
            inxst(i) = inxed(i-1)+find(stimtime(ist:ied)==1,1,'first');
        end
        rep(inxst(i):end)=i;
    end
end