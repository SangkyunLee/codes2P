function [out, inxout] = slow_xy(CMs,htwin,CM0, selframes, dir)
% function [out, inxout] = slow_xy(CMs,htwin,CM0, selframes, dir)
% CMs:centers of mass, 
% htwin: half time window
% selframes: frames to be applied
% dir: forward(1), backward(-1), bi-direction(2; default)

    if nargin<5
        dir=2;
    end

    CMx = squeeze(CMs(1,:,:));
    CMy = squeeze(CMs(2,:,:));
    if ndims(CMx) && size(CMx,1)==1,
        CMx = CMx(:);
        CMy = CMy(:);
    end

    % [~,inx] = min(abs(CMx),[],2);
    L = length(selframes);
    out = zeros(2,L);
    inxout = zeros(2,L);
    for i = 1 : L
        j = selframes(i);
        
%         if j==2108, keyboard; end
        
        
        if j> htwin
            i1 = i-htwin;
            if i1<1, i1=1; end
            nn = selframes(i1:i-1);
            inxn = j-htwin : j-1;
        else
            nn = selframes(1:i-1);
            inxn = 1 :j-1;
        end

        if j<(L-htwin)
            iL =i+htwin;
            if iL>L, iL =L; end
            np = selframes(i+1 : iL);
            inxp = j+1 : j+htwin;            
        else
            np = selframes(i+1 : L);
            inxp = j+1 : L;
        end
        if dir==1,            
            inxn = setdiff(inxn,[]);            
            inxp = [];
        elseif dir==-1 %backward
            inxn = [];
            inxp  = setdiff(inxp, []);
        elseif dir==2
            inxn = setdiff(inxn,nn);
            inxp  = setdiff(inxp, np);
        end
        

        Mn = CMx(inxn,:);
        Mp= CMx(inxp,:);

        if size(CM0,2)==1,
            CMM=CM0;
        else
            K = CM0(:,[inxn inxp]);
            K(K(:)==0)=NaN;
            CMM = nanmean(K,2);
        end
        
        mMx = mean_nCM(Mn,Mp, CMM(1));    
        if isnan(mMx)
            [~, selxi] = min(abs(CMx(j,:) -CMM(1)));
            [~, selxi2] = max(abs(CMx(j,:) -CMM(1)));
        else
            [~, selxi] = min(abs(CMx(j,:) -mMx));
            [~, selxi2] = max(abs(CMx(j,:) -mMx));
        end
        

        Mn = CMy(inxn,:);
        Mp= CMy(inxp,:);



        mMy = mean_nCM(Mn,Mp, CMM(2));
        if isnan(mMx)
            [~, selyi] = min(abs(CMy(j,:) - CM0(2)));
            [~, selyi2] = max(abs(CMy(j,:) - CM0(2)));
        else
            [~, selyi] = min(abs(CMy(j,:) - mMy));
            [~, selyi2] = max(abs(CMy(j,:) - mMy));
        end

        out(:,i)=[CMx(j,selxi), CMy(j,selyi)];
        CMx(j,selxi2)= CMx(j,selxi);
        CMy(j,selyi2)= CMy(j,selyi);
        inxout(:,i)=[selxi selyi];
    end
    CMs(1,:,:)=CMx;
    CMs(2,:,:)=CMy;
end



function mM = mean_nCM(Mn,Mp, CMk)
    if ~isempty(Mn)
        [~,inx1] = min(abs(Mn-CMk),[],2); 
        inx11 = sub2ind(size(Mn),1:size(Mn,1),inx1');    
        M = mean(Mn(inx11));
    else
        M =[];
    end
    if ~isempty(Mp)
        [~,inx2] = min(abs(Mp-CMk),[],2);
        inx22 = sub2ind(size(Mp),1:size(Mp,1),inx2');
        M = [M mean(Mp(inx22))];
    end
    if exist('M','var')
        mM = mean(M);
    else
        mM =NaN;
    end
end

