function [XY, R] = detect_lmov(XY,R,dfX, s)
% function [XY, R] = detect_lmov(XY,R,dfX, s)
% detect fast/large eyemovement such as saccades
% discard/exclude unstable pupil tracking samples
% XY, R, dfX : cell(nscan,1)--> cell(ntrial,1)
% XY: pupil center, R: pupil size in radius, dfX: paramter changes in time
% s: threshold to cut off 

    is = length(XY);
    dfX1 = cell(is,1);
    for i = 1 : is
        dfX1{i} = cell2mat(dfX{i});
    end
    dfX1 = cell2mat(dfX1);
    m0 = mean(dfX1,1);
    stdv = std(dfX1,0,1);
    
    
    
    
    for i = 1 : is
        K = dfX{i};
        XYi = XY{i};
        Ri = R{i};
        n = length(K);
        for j = 1 : n
            if isempty(K{j})
                continue;
            end
            dfXz = bsxfun(@rdivide, abs(bsxfun(@minus, K{j},m0)),stdv);
            b = dfXz>s;
            vs = ~(b(:,3)); % valid sample from distribution of pupil size changes
            
            % when too few samples, set b invalid 
            if length(vs)<0.7*size(b,1) 
                b1 = true;
            else
                b1 = any(b(vs,1) | b(vs,2)); % detect saccades
            end
            
            
            % accept only reliable and insignificant motion data
            if b1
                XYi{j}=[];
                Ri{j}=[];
            else
                XYi{j} = XYi{j}(vs,:);
                Ri{j} = Ri{j}(vs,:);
            end
             
        end
        XY{i} = XYi;
        R{i} = Ri;
    end
    
   
end