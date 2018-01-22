classdef FrameT < Timestamp
    % classdef FrameT < Timestamp
    % get_tstamp return timestamp of frame start
    
    
    
    methods
        function self = FrameT(t,y)
            self = self@Timestamp(t,y);
        end
        
        function [tinx, tstamp] = get_tstmp(self,N)
            %threshold for voltage trigger when new frame is acquired 
            % trigger (5v to 0v when new frame is acquired)
            thr = 4; 
            t = self.t;
            y = self.y;
            
            dy = [diff(y); 0];
            inx = find(dy>thr);
            % after the last frame, the voltage go back to 0V, so one more
            % count >thr are dected.
            if isempty(N)
                N= length(inx)-1;
            end
            if length(inx)<N
                error('Identified frame is %d',length(inx)-1);
            end
            tinx = [inx(1:N) inx(2:N+1)-1]; % frame start and end
            tstamp = t(inx(1:N));
            
        end
    end
end