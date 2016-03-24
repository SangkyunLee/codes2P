function self = filt(self,H)
% function self = histeq_(self)
% apply histogram equalization series of image frames
assert(isa(self,'filt_video'),'input should be filt_video');

Nframe = size(self.data,3);
for i = 1 : Nframe
    if self.KeepOrg
        self.filtdata(:,:,i) = imfilter(self.filtdata(:,:,i),H,'same');
    else
        self.data(:,:,i) = imfilter(self.data(:,:,i),H,'same');
    end
end
