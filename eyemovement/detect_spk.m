function [mdist,inx]= detect_spk(CM,htwin,thr)


len = size(CM,2);
len1 = len-htwin;
dist = zeros(htwin*2,len1);
for i = [(-htwin:-1) (1:htwin)]
    if i<0,
        j =i+htwin+1;
    else
        j =i+htwin;
    end
    
    y = CM(:,htwin+1:len1)-CM(:,(htwin+1:len1)+i);
    dist(j,htwin+1:len1)=sqrt(sum(y.^2,1));
end

mdist = mean(dist,1);
inx = find(mdist>thr);

