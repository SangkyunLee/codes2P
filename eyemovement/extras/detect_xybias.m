function [XY, R] = detect_xybias(XY,R, s)
% function [XY, R] = detect_xybias(XY,R, s)
% discard/exclude pupil samples xy-biased away from the distribution center
% XY, R: cell(nscan,1)--> cell(ntrial,1)
% XY: pupil center, R: pupil size in radius
% s: threshold to cut off 

    if nargin<2,
        s = 1.5;
    end
    is = length(XY);
    XYa = cell(is,1);
    for i = 1 : is
        XYa{i} = cell2mat(XY{i});
    end
    XYa = cell2mat(XYa);
    %m0 = mean(XYa,1);
    %stdv = std(XYa,0,1);    
    %M = bsxfun(@rdivide, abs(bsxfun(@minus, XYa,m0)),stdv);
    
    y = quantile(XYa,[.25 .75]); 

    IQ = y(2,:)-y(1,:);
    lb = y(1,:)-s*IQ;
    ub = y(2,:)+s*IQ;
    
    
    for i = 1 : is
    
        XYi = XY{i};
        Ri = R{i};
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
        XY{i} = XYi;
        R{i} = Ri;
    end
    
   
end