function self = threshold(self,thr,op)
% function self = threshold(self,thr,op)
% when op='up', data(data(:)<thr)=thr;
% when op='down', data(data(:)>thr)=thr;
assert(isa(self,'filt_video'),'input should be filt_video');

Nframe = size(self.data,3);
for i = 1 : Nframe
    if self.KeepOrg
        tmpImg = self.filtdata(:,:,i);                
    else
        tmpImg = self.data(:,:,i);
    end
    
    switch op
        case {'down'}
            tmpImg(tmpImg(:)<thr)=thr;
        case{'up'}
            tmpImg(tmpImg(:)>thr)=thr;
    end                
    
    if self.KeepOrg
        self.filtdata(:,:,i) = tmpImg;
    else
        self.data(:,:,i) = tmpImg;
    end
    
end
