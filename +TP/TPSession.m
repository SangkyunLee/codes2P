classdef TPSession 
    % 2016-02-06, Sangkyun Lee
    
    properties
        fn_data; % data_file after calling load_data
        
       
        session_path; % root_path for all session data
        
        animal_state; % 'AW','AN(fentanyl)', etc
        imagingmethod; % spiral(224), sprial(512), resonant(1024)
        
        no_scan; % # scans
        scanId; % scan identity numbers
        scan_stim; % stimulus type in each scan        
        scan_ncell; % no cell identified
        FOV_size; % FOV size in pixel   
        
        
        
        scans=struct([]); % data loaded by load_data function
        motion=struct([]); % motion_estimate (unit: image frame)        
        eye = struct([]); 
        
    end
    methods 
        function self =TPSession(fn_data,  info)
            if exist(fn_data,'file')~=2
               error('File not exist:%s',self.fn_data);
               % TO DO for later % loading 
               % I have to make a function to load data into a matlabfile
               % load_data(info)
            else
                loaddata =load(fn_data);   
                
                
                self.scans = loaddata.data;                
                self.fn_data = fn_data;
                self.no_scan = length(loaddata.data);  
                self.scanId = zeros(1,self.no_scan);
                for iscan= 1: self.no_scan
                    snum = strsplit(loaddata.data(iscan).Params.files.subpath_xml,'-');
                    self.scanId(iscan) =  str2double(snum{end});
                end
                
                self.scan_ncell = zeros(self.no_scan,1);
                if nargin>1,     
                    
                    pnames = fieldnames(info);
                    for ip = 1 : length(pnames)
                        self.(pnames{ip})=info.(pnames{ip});
                    end
                    
                end
                self.FOV_size = size(self.scans(1).meanimg);
                
                for iscan = 1: self.no_scan
                    self.scan_stim{iscan}= self.scans(iscan).exptype;
                    ncell = size(self.scans(iscan).Nhat,2);
%                     if length(self.scans(iscan).ROI) > 1,
%                         for iseg = 2 : length(self.scans(iscan).ROI)
%                             ncell = min(ncell,length(self.scans(iscan).ROI{iseg}));
%                         end
%                     end
                    self.scan_ncell(iscan) = ncell;
                end
                clear loaddata;               
            end            
        end
        
        
        
        function data = get_data(self)
            data=self.scans;
        end
        
        
        function [b bs] = check_stim_consistency(self, scanlist,evtsort)
            % function [b bs] = check_stim_consistency(self, scanlist,evtsort)
            % if consistent, return true, reference point is the
            % scanlist(1)
            
            firstscan = scanlist(1);
            stimparam1 = self.scans(firstscan).Params.stimparam;
            
            bs = ones(length(scanlist),length(evtsort));
            for ievt = 1 : length(evtsort)                
                evtlist = stimparam1.(lower(evtsort{ievt}));
                for inx_scan0 = 1:length(scanlist)-1
                    inx_scan = scanlist(inx_scan0);
                    stimparam = self.scans(inx_scan).Params.stimparam;                    
                    evtlist2 = stimparam.(lower(evtsort{ievt}));
                    if length(evtlist) == length(evtlist2)
                        bs(inx_scan0,ievt) = ~ sum(abs(evtlist-evtlist2));
                    else
                        bs(inx_scan0,ievt)= false;
                        
                    end
                end
                
            end
     
            if all(bs(:)),
                b=true;
            else
                b=false;
            end
        end
        function [b bs] = check_param_consistency(self, scanlist,parsort)
            % bs = check_param_consistency(self, scanlist,parsort)
            % if consistent, return true, reference point is the
            % scanlist(1)
            
            firstscan = scanlist(1);
            param1 = self.scans(firstscan).Params;
            
            bs = ones(length(scanlist),length(parsort));
            for ipar = 1 : length(parsort)                
                parlist = param1.(parsort{ipar});
                if strcmp(parsort{ipar},'msperframe')
                    parlist = floor(parlist);
                end
                    
                for inx_scan0 = 1 : length(scanlist)-1%scanlist(2:end)
                    inx_scan = scanlist(inx_scan0);
                    param = self.scans(inx_scan).Params;                    
                    parlist2 = param.(parsort{ipar});
                    if strcmp(parsort{ipar},'msperframe')
                        parlist2 = floor(parlist2);
                    end 
                    
                    if length(parlist) == length(parlist2)
                        bs(inx_scan0,ipar) = ~ sum(abs(parlist-parlist2));
                    else
                        bs(inx_scan0,ipar)= false;
                        
                    end
                end
                
            end
            
            if all(bs(:)),
                b=true;
            else
                b=false;
            end
        end
        
        
        function self = set_motion_rotaroad_imageframe(self,twin,channel)
            % function set_motion_rotaroad_imageframe(self,twin)
            % calculate motion parameters and set the motion paramters
            for iscan = 1 : self.no_scan
                files = self.scans(iscan).Params.files;
                DAQpath = fullfile(self.session_path,files.subpath_xml);
                FMOTpath = DAQpath;
                frame_period = self.scans(iscan).Params.msperframe/1000;
                pixel_resolution = self.scans(iscan).ROI{1}(1).pixelres;
                params=struct('DAQpath',DAQpath,...
                    'FMOTpath',FMOTpath,...
                    'frame_period',frame_period,...
                    'pixel_resolution',pixel_resolution,...
                    'twin',twin,...
                    'channel',channel);
                [~, dist_twin, IMG_motion, ~]= TP.TPSession.extract_motion_combo1(params);
                
                % frame_start in DAQ signal
                frame_start = self.scans(iscan).timeinfo.frame_start;                
                dist_twin_frame = dist_twin(frame_start(1:length(dist_twin))>0);
                if length(dist_twin_frame)<length(IMG_motion)
                    dist_twin_frame1 = ones(size(IMG_motion))*500;
                    dist_twin_frame1(1:length(dist_twin_frame))= dist_twin_frame;
                    dist_twin_frame = dist_twin_frame1;
                else
                    dist_twin_frame = dist_twin_frame(1:length(IMG_motion));
                end
                
                motion1.IMG_motion = IMG_motion;                
                motion1.dist_twin_frame = dist_twin_frame; 
                if isempty(self.motion)
                    self.motion =motion1;                    
                else
                    self.motion(iscan) =motion1;
                end
            end
        end
        
        function self = set_eyepar(self,par)
            % function self = set_eyepar(self,par)
            % set pupil parameters(x,y, size in radius) in trials
            
 
            ns = self.no_scan;
            eye0 = struct('PP',[],'xy',[],'r',[],'vf_trial',[]);
            self.eye = repmat(eye0,[1 ns]);
            

            dfX = cell(self.no_scan,1);
            for i = 1 : ns
                Params = self.scans(i).Params;
                timeinfo = self.scans(i).timeinfo;
                Params.files.mainpath = self.session_path;
                [XYi, Ri, dPi, PP, vfinxi] = TP.get_eyepar_trial(Params,timeinfo);
                self.eye(i).xy = XYi;    
                self.eye(i).r = Ri;    
                self.eye(i).PP = PP;
                self.eye(i).vf_trial = vfinxi;
                dfX{i} = dPi;                
            end
            
%             self.eye = TP.detect_lmov(self.eye,dfX, par.lmovthr); % detect the large eye movement (e.g., saccades)
%             self.eye = TP.detect_xybias(self.eye,par.dxythr); % detect the pupil center away from the pupil center distr           
%             self.eye = TP.detect_slpup(self.eye,par.puppar); % detect small and large pupil
            
            eye1 = TP.detect_lmov(self.eye,dfX, par.lmovthr); % detect the large eye movement (e.g., saccades)            
            eye2 = TP.detect_xybias(self.eye,par.dxythr); % detect the pupil center away from the pupil center distr                       
            eye3 = TP.detect_slpup(self.eye,par.puppar); % detect small and large pupil
            ntr = zeros(ns,1);
            for i = 1: ns
                r1 = cellfun(@isempty,eye1(i).r);
                r2 = cellfun(@isempty,eye2(i).r);
                r3 = cellfun(@isempty,eye3(i).r);
                iival = find(r1 | r2 | r3);
                emptycell = cell(length(iival),1);
                eye1(i).xy(iival)=emptycell;
                eye1(i).r(iival)=emptycell;
                eye1(i).vf_trial(iival)=emptycell;            
                ntr(i)=240-length(iival);
            end
            fprintf('Total trial= %d: %s trials selected in scans\n',sum(ntr),num2str(ntr')); 
            self.eye = eye1;
        end
        
        
        
        function [X, Xc, seltrial, timeF] = collect_trialresp_singlescan(self, inx_scan, pre_stimtime_ms, motionthr, beye)
            % function [X, Xc] = collect_trialresp_singlescan(self, inx_scan, pre_stimtime_ms, motionthr, beye)
            % This function works only for trial-based experiment
            
            if nargin<4,
                motionthr = struct([]);
            end
            if nargin<5
                beye = false;
            end
            
            param = self.scans(inx_scan).Params;
            nframe_prestim = round(pre_stimtime_ms/param.msperframe);
            
            stimparam = param.stimparam;
            
            stimtime = stimparam.stim_samplesinNI/param.samplingfreq_NI;
            blanktime = stimparam.blank_samplesinNI/param.samplingfreq_NI;

            trialtime = stimtime+blanktime;
            nframe_trial = round(trialtime/param.msperframe*1000);            
            
            timeF = (-nframe_prestim:nframe_trial-1)*param.msperframe; % in millisecond
            
            
            spec.frames = -nframe_prestim : nframe_trial-1;
            spec.dataType='dFF';
            spec.nCell= self.scan_ncell(inx_scan);
            
            if isstruct(motionthr) && ~isempty(motionthr)
                motion1.tmotion1 = self.motion(inx_scan).IMG_motion;
                motion1.tmotion2 = self.motion(inx_scan).dist_twin_frame;
                motion1.motionthr1 = motionthr.IMG_motion_thr;
                motion1.motionthr2 = motionthr.dist_twin_thr;
                motion1.op='|';
                motion1.bdisp = false;
                spec.motion = motion1;
            end
            Xc = data_sort(self.scans(inx_scan),spec,1);
            
            spec.dataType='Nhat';    
            [X, others]=data_sort(self.scans(inx_scan),spec,1);
            
            seltrial = others.seltrial{1};
            if beye
                eyei = self.eye(inx_scan);            
                eyetrial = find(~cellfun(@isempty,eyei.xy));                
                inxmap = zeros(max(max(seltrial),max(eyetrial)),1);
                inxmap(seltrial)=1:length(seltrial);                
                seltrial1 = intersect(seltrial,eyetrial);
                
                newinx = inxmap(seltrial1);
                Xc{1} = Xc{1}(newinx,:,:);
                X{1} = X{1}(newinx,:,:);
               
                seltrial = seltrial1;
            end
        end
        
        function CL = get_cell_spkestval(self, scanlist)
            % function CL = get_cell_spkestval(self, scanlist)
            % return mask of cell for successful spike estimation in
            % scanlist
            ncell = self.scan_ncell(1);
            CL =ones(1,ncell);
            for inx_scan = scanlist
                if isfield(self.scans(inx_scan),'info') && isfield(self.scans(inx_scan).info,'history')
                    for icell = 1:ncell
                        CL(icell) = CL(icell) &self.scans(inx_scan).info.history(icell);
                    end
                else
                    CL1 = sign(sum(self.scans(inx_scan).Nhat,1));
                    for icell = 1:ncell
                        CL(icell) = CL(icell) & CL1(icell);
                    end
                end  
            end
            
        end
        
        
        function [MXspk, MXc, trialinxs, tf] = collect_trialMresp(self, scanlist, twins ,motionthr, beye)
            % function collect_trialMresp(self, scanlist, twins ,motionthr,beye)
            % twins:    time windows to calculate the mean response
            %           can be multiple segments with matrix 2x no. window
            %           twins(1,i): start time (millisecond) in window i
            %           twins(2,i): end time (millisecond) in window i
            %           0 is the reference time synchronized with stimulus
            %           onset
            
            if nargin<4,
                motionthr = struct([]);
            end
            if nargin<5,
                beye = false;
            end
            [d nwin] = size(twins);
            assert(d==2,'twins should be a matrix with a size of 2xnwindow');
            assert(max(scanlist)<= self.no_scan, ['Total number of scan is ' num2str(self.no_scan)]);
            
            
            pre_stimtime_ms = min(twins(1,:));
            if pre_stimtime_ms >0
                pre_stimtime_ms = 0;
            else
                pre_stimtime_ms = -pre_stimtime_ms; % should be positive for collect_trialresp_singlescan
            end
            
            ns = length(scanlist);
            MXspk = cell(ns,nwin);
            MXc = cell(ns,nwin);
            tfs = cell(ns,1);
            trialinxs = cell(ns,1);
            
            
            for inx_scan0 = 1:ns
                %param = self.scans(inx_scan).Params;
                %stimparam = param.stimparam;
                inx_scan = scanlist(inx_scan0);
               
                [X, Xc, seltrial, tf] = collect_trialresp_singlescan(self, inx_scan, pre_stimtime_ms, motionthr, beye);
               
                
                
                for iwin = 1: nwin
                    inxframes = (tf>twins(1,iwin) & tf<=twins(2,iwin));  
                                       
                    mXspk = squeeze(mean(X{1}(:,inxframes,:),2));                        
                    mXc = squeeze(mean(Xc{1}(:,inxframes,:),2));                                            
                    MXspk{inx_scan0,iwin} = mXspk;
                    MXc{inx_scan0,iwin} = mXc;
                end
                
                trialinxs{inx_scan0} = seltrial;
                tfs{inx_scan0} = tf;
                
            end
            
        end
        
        
        
        
       
    end
    
    
    
    
    
    
    
    methods (Static )
        % 
        function [speed, dist_twin, IMG_motion, outparam]= extract_motion_combo1(params)
            %function [speed, dist_twin, IMG_motion, outparam]= extract_motion_combo1(params)
            % extract motion parameters from rotaroad and 2P image
            
            
            DAQpath = params.DAQpath;  % DAQ for rotaroad rotation
            FMOTpath = params.FMOTpath; % image motion parameter path
            pixel_resolution = params.pixel_resolution;
            frame_period = params.frame_period; % TP image frame period (
            twin = params.twin; % time window for movement distance 
            channel = params.channel; % DAQ channel index
            
            [speed, dist_twin, IMG_motion, outparam]  =...
                extract_motion_combo(DAQpath,FMOTpath,pixel_resolution, frame_period,twin,channel);
        end
    end
end