function self = histeq_(self)
% function self = histeq_(self)
% apply histogram equalization series of image frames
assert(isa(self,'filt_video'),'input should be filt_video');

Nframe = size(self.data,3);
for i = 1 : Nframe
    if self.KeepOrg
        self.filtdata(:,:,i) = histeq(self.filtdata(:,:,i));
    else
        self.data(:,:,i) = histeq(self.data(:,:,i));
    end
end
