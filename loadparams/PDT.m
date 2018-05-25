classdef PDT < Timestamp
    % classdef PDT < Timestamp
    % get_tstamp return timestamp of photodiode_event

    
    methods
        function self = PDT(t,y)            
            self = self@Timestamp(t,y);
        end
        
        function [tinx, tstamp] = get_tstmp(self,minint,mode)
            % function [tinx, tstamp] = get_tstmp(self,minint,mode)
            % minint: min interval between two photodiode-event
            % mode: 'auto' vs 'manual'
            
            t = self.t;
            y = self.y;                        
            tinx = get_evt(t,y,minint,mode);
            tstamp = t(tinx);         
            
        end
        
        function plot_result(self, tstamp, t1, title)
            
            if nargin<3
                t1 = 40;
            end
            if nargin<4
                figure;
            else
                figure('Name',title);
            end
            subplot(2,1,1);
            inx1 = self.t<t1;
            inx2 =find(tstamp<t1);
            plot(self.t(inx1),self.y(inx1)); 
            hold on; plot(tstamp(inx2),max(self.y)*ones(size(inx2)),'r.');
            subplot(2,1,2);
            inx1 = self.t>self.t(end)-t1;
            inx2 =find(tstamp>(self.t(end)-t1));
            plot(self.t(inx1),self.y(inx1)); 
            hold on; plot(tstamp(inx2),max(self.y)*ones(size(inx2)),'r.');
        end
    end
end

function inxt = get_evt(t,y,minint,mode)
    if strcmp(mode,'auto')
        thr = mean(y)+5*std(y,0); 
        inxt0 = 1;        
    else
        len = round(length(y)*0.2);
        h=figure; plot(t(1:len,1),y(1:len)); xlabel('time (sec)')
        waitfor(h)
        
        prompt = {'Threshold','Time start'};
        dlg_title = 'Input a starting point for the photo-diode';
        num_lines = 2;
        def = {'3','5'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        thr = eval(answer{1});
        t0 = eval(answer{2});
        inxt0 = round(t0*1/diff(t(1:2)))+1;
        
        t = t(inxt0:end);
        y = y(inxt0:end);        
    end    
    inx = find(y>thr);
    dt = diff(t(inx));     
    inx2 =[1; find(dt> 0.7*minint)+1];
    inxt = inx(inx2);      
    inxt = inxt+inxt0-1;

    
end
