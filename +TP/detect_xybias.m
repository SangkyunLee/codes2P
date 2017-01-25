function eye = detect_xybias(eye, s)
% function eye = detect_xybias(eye, s)
% discard/exclude pupil samples xy-biased away from the distribution center
% eye().xy, eye().r:  cell(ntrial,1)
% xy: pupil center, r: pupil size in radius
% s: threshold to cut off 

    if nargin<2,
        s = 1.5;
    end
    is = length(eye);
    XYa = cell(is,1);
    for i = 1 : is
        % generate artifical samples for invalid trials
        xy = eye(i).xy;
        samplelen = cellfun(@length,xy);        
        iival = find(samplelen==0);
        if ~isempty(iival) && any(samplelen>0)
            mlen = round(mean(samplelen(samplelen>0)));
            tmp={NaN*ones(mlen,2)};
            tmp = repmat(tmp,[length(iival),1]);
            xy(iival)= tmp(:);
        end
        %----------------------
        
        XYa{i} = cell2mat(xy);
    end
    XYa = cell2mat(XYa);
    
    % fill outliers in invalid trials
    [lb0, ub0] = get_bound(XYa,10);  
    inx_nan = find(isnan(XYa(:,1)));
    L= length(inx_nan);
    i2f = floor(L/2);
    i2s = L-i2f;
    k = randperm(L);
    XYa(inx_nan(k(1:i2f)),:)= ones(i2f,1)*lb0;
    XYa(inx_nan(k(i2f+1:end)),:)= ones(i2s,1)*ub0;
    %------------------------------------

    
    [lb, ub] = get_bound(XYa,s);
    
    
    for i = 1 : is
    
        XYi = eye(i).xy;
        Ri = eye(i).r;
        n = length(XYi);
        for j = 1 : n
            if isempty(Ri{j})
                continue;
            end
            
            
            % detect outlier
            b = bsxfun(@minus,XYi{j},lb)<0 | bsxfun(@minus,XYi{j},ub)>0;
            b = b(:,1) | b(:,2);
            
            % true if close to distribution median
            b = ~b;
            
            % set invalid when less than 50% samples are away from the sample median         
            if sum(b)< 0.5*length(b)
                XYi{j}=[];
                Ri{j}=[];
            else
                XYi{j} = XYi{j};
                Ri{j} = Ri{j};
            end
             
        end
        eye(i).xy = XYi;
        eye(i).r = Ri;
    end
    
   
end

function [lb, ub] = get_bound(X,s)
    y = quantile(X,[.25 .75]); 

    IQ = y(2,:)-y(1,:);
    lb = y(1,:)-s*IQ;
    ub = y(2,:)+s*IQ;
end