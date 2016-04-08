function [out, inx_miss, inx_fill] = find_appCM(In,winsize,minN)

CM =In.CM;
sPIX = In.sPIX;
out = In;
inx_miss = find(sum(CM,1)==0);
inx_fill = zeros(size(CM,2),1);
hwin = floor(winsize/2);
for i=inx_miss
    wini = i-hwin:i+hwin;
    if wini(1)<1 ||wini(end)>size(CM,2)
        continue;
    end
    CMi = CM(:,wini);
    inx2 = find(sum(CMi,1)~=0);
    if length(inx2)>minN,
        out.CM(:,i)=mean(CMi(:,inx2),2);
        for j = inx2
            if j== inx2(1)
                out.sPIX{i} = sPIX{wini(j)};
            else
                out.sPIX{i} = intersect(out.sPIX{i}, sPIX{wini(j)});
            end
        end
        inx_fill(i) = 1;
    end
   
end
inx_fill = find(inx_fill==1);
inx_miss = setdiff(inx_miss, inx_fill);
        