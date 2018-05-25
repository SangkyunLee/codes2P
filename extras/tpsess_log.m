classdef tpsess_log
    properties
        info;
        nscan;
    end
    
    methods
        function self = tpsess_log(logfn,sheet)
            [~,~,raw] = xlsread(logfn,sheet);

            finfo=struct([]);
            nfld = size(raw,2);
            flds = cell(1,nfld);
             irow=1;
            for icol=1:size(raw,2)
                if ~isnan(raw{irow,icol})            
                    flds{icol} = raw{irow,icol};       
                end
            end

            kk=1;
            for irow=2:size(raw,1)    
                if  ~isnan(raw{irow,1})
                    for icol=1:length(flds)
                        if ~isempty(flds{icol})
                            finfo(kk).(flds{icol}) = raw{irow,icol};            
                        end
                    end
                    kk = kk+1;
                end
            end
            
            
            self.info = finfo;
            self.nscan = length(finfo);
        end
        
        function list = search_scan(self,varargin)
           
            
            val = true*ones(self.nscan,1);
            
            for i = 1 : self.nscan
                val(i) = self.check(i,varargin);
            end
            list = find(val);            
        end
    end
      
    
    
    
    methods (Access=private)        
        function val = check(self, iscan,varargin)
            finfo = self.info(iscan);
            if iscell(varargin) && length(varargin)==1
                varargin = varargin{1};
            end
            
            nfld = floor(length(varargin)/2);            
            for i = 1: nfld
                if i==1,
                    val =  strcmp(finfo.(varargin{2*(i-1)+1}),varargin{2*i});
                else
                    val = val & strcmp(finfo.(varargin{2*(i-1)+1}),varargin{2*i});
                end
            end
            
            
        end
    end
end

