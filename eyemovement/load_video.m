function self = load_video(videofn,frame,Xr,Yr)

% self =read_subsample(self, [],Xr,Yr);

self= filt_video(videofn);
if nargin<3
    self =read_subsample(self, frame);
else
    self =read_subsample(self, frame, Xr, Yr);
end
heq = @(self) histeq_(self);
H=fspecial('gaussian',5,5);
sF = @(self) filt(self,H);
self = addfiltstep(self,heq);
self = addfiltstep(self,sF);
applyfilt(self); 