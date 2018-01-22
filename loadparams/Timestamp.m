classdef Timestamp
    
    properties 
        t; % time        
        y; % signal

        
    end
    methods
        function self = Timestamp(t,y)
            assert(length(t)==length(y));
            self.t = t;
            self.y =y;
        end
        
        function self = get_tstmp(self)
        end
        
        
        function v = get(self,fld)
            v = self.(fld);           
        end
    end
end