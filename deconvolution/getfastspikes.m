function [Nhat, As, Betas,historys]= getfastspikes(F,sT, opts)
% function [Nhat trends]= getfastspikes(F,sT, opts)
% INPUT:
%     F: TxC raw fluorescence matrix or 1xC fluorescence cell
%     sT: chunck length
%     opts.dt (frame duration in sec)
%     opts.tau 1(default)
%     opts.ndct: drift correction
%     opts.H: regressor matrix
%     opts.F0_thr: F0 threshold
% OUTPUT:
%     Nhat
%     Beta: if ndct=0, zero trends
%
% 2013-10-30 Sangkyun Lee
% 2013-12-04 modifiedy by Sangkyun Lee
%     add deconvolution in multiple pixels within a cell

if iscell(F)
    nCells = length(F);
    [T11 T12] = size(F{1});
    [T21 T22] = size(F{2});
    if T11 == T21, 
        T=T11;
    elseif T12 == T22
        T=T12;
    else
        error('mismatch in no. frame');
    end
else
    [T nCells] = size(F);
end
Nhat=zeros(T,nCells);
if isempty(sT)
    nchunck =1;
    sT=T;
else
    nchunck = round(T/sT);
end
As = cell(nchunck,nCells);
Betas = cell(1,nCells);
historys = cell(nchunck,nCells);
if isfield(opts,'ndct')
    ndct = opts.ndct;
    tc = linspace(0,2*pi,T)';    
    ndct1 = 2*ndct+1;
    dct(1:T,1:ndct1) = cos(tc*(ndct:-0.5:0));
else
    ndct=0;
    dct =[];
end

if isfield(opts,'F0_thr')
    F0_thr=opts.F0_thr;
else
    F0_thr =0.2;
end

if isfield(opts,'H') && ~isempty(opts.H)
    H1 = opts.H;
else
    H1=[];
end

if ndct>0
    H2 = dct;
else
    H2 = [];
end

H = [H1 H2];


for ic=1:nCells
    if iscell(F)
        Fi = F{ic};
        if size(Fi,2)==T
            Fi = Fi';
        end
    else        
        Fi = F(:,ic);    
    end
    
    if ~isempty(H)
        Betas{ic} = H\Fi;    
        Fi = Fi - H(:,1:end-1)*Betas{ic}(1:end-1,:);  
    end

    v = sort(mean(Fi,2));
    F0 = mean(mean(v(1:round(T*F0_thr),:)));
    Y0 = (Fi -F0)/F0;
%     Y0 = bsxfun(@rdivide,bsxfun(@minus,Fi,F0),F0);
    
    
    
    for ichunck=1:nchunck
        tinx1 = 1+(ichunck-1)*sT;
        if ichunck>1,
            tinx1 = tinx1 -100;
        end
        if ichunck == nchunck
            tinx2 = T;
        else
            tinx2 = ichunck*sT;
        end
        Y =Y0(tinx1:tinx2,:);       
        [npart,a,history] = getspikes(Y,opts);
       
        if ichunck==1,
            Nhat(tinx1:tinx2,ic) = npart(1:sT);        
        else
            Nhat(tinx1+100:tinx2,ic) = npart(101:end);            
        end 
        As{ichunck,ic}=a;
%         history = rmfield(history,'ns');
%         history = rmfield(history,'as');        
        historys{ichunck,ic}=history;
    end
end
