classdef filt_video < VideoReader

    properties
        data; % video frames loaded, this is only for acrhomatic movie
        filtdata; % filtered data
        frames; % frames loaded
        Xr; % x-range
        Yr; % y-range
        filtsteps; % filtersteps
        KeepOrg; % boolean to save original data
        bfilted; % boolean to indicate whether filters were applied into data
    end

    methods
        function self = filt_video(fnvideo, KeepOrg)
            self = self@VideoReader(fnvideo);
            self.bfilted = false;
            if nargin<2
                self.KeepOrg = false;
            else
                self.KeepOrg = KeepOrg;
            end
        end

        
        function self = read_subsample(self, frames, Xr, Yr)
        % function self = read_subsample(self, frames, Xr, Yr)
        % loading video segments
            if nargin<3
                Xr(1)=1; Xr(2)=self.Width;
                Yr(1)=1; Yr(2)=self.Height;    
            end
            if Xr(1)<1, Xr(1)=1; end
            if Xr(2)>self.Width, Xr(2)=self.Width; end
            if Yr(1)<1, Yr(1)=1; end
            if Yr(2)>self.Height, Yr(2)=self.Height; end
            self.Xr = Xr;
            self.Yr = Yr;



            if isempty(frames)
                frames = 1:self.NumberOfFrames;
                self.frames = frames;                
            else
                self.frames = frames;                
            end
            Nframe = length(frames);
                
            for ifr = 1 : Nframe
                frame  = frames(ifr);
                if frame == frames(1),
                    dtype = check_videotype(self);                          
                    eval(sprintf('self.data = %s(zeros(Yr(2)-Yr(1)+1,Xr(2)-Xr(1)+1,Nframe))',dtype));                    
                end
                tmpImg = read(self, frame);
                self.data(:,:,ifr) = tmpImg(Yr(1):Yr(2),Xr(1):Xr(2));                            

            end
            if self.KeepOrg
                self.filtdata = self.data;
            end
        end


        
        function dtype = check_videotype(self)
        %----- check videoformat
            switch(self.VideoFormat)
                case {'Mono8'}
                    dtype = 'uint8';
                case {'Mono12','Mono16'}
                    dtype ='uint16';
                otherwise
                    error('it is not implemented yet');
            end
        end

        function self = addfiltstep(self,step)
        %-----add filtersteps
            index = numel(self.filtsteps)+1;
            self.filtsteps{index}=step;
        end
  
        function self = applyfilt(self)   
        %----- apply filters
            for i = 1 : numel(self.filtsteps)
                self = self.filtsteps{i}(self);
            end
            self.bfilted = true;
                
        end

    end

end