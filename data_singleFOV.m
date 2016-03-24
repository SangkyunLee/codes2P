classdef data_singleFOV
    properties (SetAccess = private)
        fn_data; % datafile name; the data specified by this file name should be obtained by "loaddata.m"
        data; % data structure saved in fn_data; structure array, whose component consists of a single scan data        
        bawake = false;
        motionfn={}; % when awake
        

        
%         MRs=struct('mC',[],'mN',[],'ori',[],'contrast',[],'trial_ori',[],...
%             'trial_cont',[],'avgresp_spk',[],'avgresp_dFFc',[],'timeTR',[],'events',[]);
        
%         MR.mC = mXc;
%         MR.mN = mXspk;
%         MR.ori = ORI_deg;
%         MR.contrast = Cont;
%         MR.trial_ori = trial_ori;
%         MR.trial_cont = trial_cont;
%         MR.avgresp_spk = avgresp_spk;
%         MR.avgresp_dFFc = avgresp_dFFc;
%         MR.timeTR = tf;
%         MR.events = events;
       
        description='';
        Nscan=0;        
    end
    
%     properties (SetAccess = private, Dependent = true)
%         
%     end
    
    methods
        function self = data_singleFOV(fn_data,motionfn, spkthr,bnorm)
            if exist(fn_data,'file')~=2
               error('File not exist:%s',self.fn_data);
            else
                loaddata =load(fn_data);   
                self.data = loaddata.data;                
                self.fn_data = fn_data;
                self.Nscan = length(loaddata.data);
                if nargin>2
                    for iscan = 1 : self.Nscan
                        Nhat = self.data(iscan).Nhat;
                       
                        maxv1 = spkthr*std(Nhat,[],1);                          
                        inxvalidcell = find(maxv1>0);
                        if nargin>3 && bnorm,                                                                                                   
                            m1 = bsxfun(@minus,Nhat, maxv1);
                            Npbin = sum((m1>0),1);
                            normvalue = zeros(1,size(m1,2));
                            for icell=1:size(m1,2)
                                normvalue(icell) = sum(m1(m1(:,icell)>0,icell))/Npbin(icell);
                            end
                            
                            Nhat(:,inxvalidcell) = bsxfun(@rdivide, Nhat(:,inxvalidcell), normvalue(:,inxvalidcell));
                            
                        end 
                        NhatP = Nhat(:,inxvalidcell);
                        stdNhat = std(NhatP(:),0,1);
                        btmp = sign(bsxfun(@minus,Nhat,spkthr*stdNhat));
                        Nhat(btmp(:)<0)=0;
                        %self.data(iscan).Nhat= log(1+Nhat*100);
                        self.data(iscan).Nhat= Nhat;
                    end
                end
                
                clear loaddata;
                if iscell(motionfn)
                    if isstr(motionfn{1}) && length(motionfn)>self.Nscan
                        self.bawake = true; 
                        self.motionfn = motionfn;
                    else
                        self.bawake = false; 
                    end
                else                                       
                    self.bawake = false; 
                end
%                 nMRs = MRs;
%                 for iscan=1:self.Nscan
%                     nMRs(iscan)=MRs;
%                 end
%                 self.MRs=nMRs;
            end            
        end
        
        
        
        function MR = get_meanresp_insinglescan(self,inxscan,motionspec)
            %  MR = collect_meanresp_insinglescan(self,inxscan)
            % collect mean response in dFF and estimated spike in a single movie
            
            if inxscan>self.Nscan
                error('Total scan: %d',self.Nscan);
            end
            
            bquiet = true;                                
            if isfield(self.data(inxscan).Params.stimparam,'stim1_samplesinNI')
                stimtime = self.data(inxscan).Params.stimparam.stim1_samplesinNI/self.data(inxscan).Params.samplingfreq_NI;
                blanktime = self.data(inxscan).Params.stimparam.blank1_samplesinNI/self.data(inxscan).Params.samplingfreq_NI;
            else
                stimtime = self.data(inxscan).Params.stimparam.stim_samplesinNI/self.data(inxscan).Params.samplingfreq_NI;
                blanktime = self.data(inxscan).Params.stimparam.blank_samplesinNI/self.data(inxscan).Params.samplingfreq_NI;
            end
            trialtime = stimtime+blanktime;
            nframe_trial = round(trialtime/self.data(inxscan).Params.msperframe*1000);
            nCell=size(self.data(inxscan).dFF,2);
            spec.frames = (-3:nframe_trial-1);
            spec.dataType = 'dFF';
            spec.nCell = nCell;
            if self.bawake
                if nargin>2 
                    load(self.motionfn{inxscan});
                    tmotionpix = sqrt(M(3,:).^2 + M(4,:).^2);
                    spec.motion.tmotion1 = tmotionpix;
                    spec.motion.motionthr1 = motionspec.motionthr1;
                    spec.motion.tmotion2 = [0 abs(diff(tmotionpix))];
                    spec.motion.motionthr2 = motionspec.motionthr2;
                    spec.motion.op='|';
                    spec.motion.bdisp = false;          
                end
            end

            Xc = data_sort(self.data(inxscan),spec,bquiet);
            clear spec;

            spec.frames=(-3:nframe_trial-1);
            spec.dataType='Nhat';
            spec.nCell=nCell;
            if self.bawake
                if nargin>2 
                    load(self.motionfn{inxscan});
                    tmotionpix = sqrt(M(3,:).^2 + M(4,:).^2);
                    spec.motion.tmotion1 = tmotionpix;
                    spec.motion.motionthr1 = motionspec.motionthr1;
                    spec.motion.tmotion2 = [0 abs(diff(tmotionpix))];
                    spec.motion.motionthr2 = motionspec.motionthr2;
                    spec.motion.op='|';
                    spec.motion.bdisp = false;          
                end
            end
            [Xspk others]=data_sort(self.data(inxscan),spec,bquiet);
            events = others.events{1};

            

            mXc1 = squeeze(mean(Xc{1},1));                
            smXc1 = std(mXc1,0,1);

            % cell selection procedure
            % select cell activity >0 in df/f and std(meanresp)< 3sigma
            learningvalid =zeros(1,nCell);
            for icell = 1:nCell
                learningvalid(icell) = ~isnan(self.data(inxscan).info.nhat_history(icell).posts(end));
            end
            inxincell = find(mean(mXc1,1)>0& smXc1<(mean(smXc1)+3*std(smXc1)) & learningvalid);
            inxoutcell = setdiff((1:nCell),inxincell);


            Xc1 = permute(Xc{1},[3 2 1]);   
            Xspk1 = permute(Xspk{1},[3 2 1]);                   

            avgresp_spk = mean(mean(Xspk1(inxincell,:,:),3),1);            
            avgresp_dFFc = mean(mean(Xc1(inxincell,:,:),3),1);
            tf = (-3:nframe_trial-1)*self.data(inxscan).Params.msperframe/1000;
            %figure; plot(tf,avgresp_dFFc,'.-');

            %- peakframes from estimated spikes
            %[~, mi1]=max(avgresp_spk);
            %peakframes_spk = mi1-1:mi1+1;
            peakframes_spk = find(tf>=0 & tf<=stimtime);


            [~, mi2]=max(avgresp_dFFc);
            %peakframes_dFFc=mi2-1:mi2+1;
            ntf =tf-tf(mi2);
            peakframes_dFFc = find(abs(ntf)<stimtime);

            %-----------------------------------------

                        
            if isfield(self.data(inxscan).Params.stimparam,'orientations')
                ORI_deg = self.data(inxscan).Params.stimparam.orientations;
            else
                ORI_deg = self.data(inxscan).Params.stimparam.orientation;                    
            end
            if isfield(self.data(inxscan).Params.stimparam,'contrasts')
                Cont = self.data(inxscan).Params.stimparam.contrasts;
            else
                Cont = self.data(inxscan).Params.stimparam.contrast;
            end

            mXspk = squeeze(mean(Xspk{1}(:,peakframes_spk,:),2));    
            mXc = squeeze(mean(Xc{1}(:,peakframes_dFFc,:),2));    
            mXspk(:,inxoutcell)=0;
            mXc(:,inxoutcell)=0;
            trial_ori = self.data(inxscan).Params.stimparam.cond.StimEventLUT(:,3);
            trial_cont = self.data(inxscan).Params.stimparam.cond.StimEventLUT(:,1);
            
            
            MR.mC = mXc;
            MR.mN = mXspk;
            MR.ori = ORI_deg;
            MR.contrast = Cont;
%             MR.trial_ori = trial_ori;
%             MR.trial_cont = trial_cont;
            MR.avgresp_spk = avgresp_spk;
            MR.avgresp_dFFc = avgresp_dFFc;
            MR.timeTR = tf;
            MR.events = events;
            
        end
        
        
        %% function self = get_meanresp_acrossscans(self, inxscans)
        function MRs = get_meanresp_acrossscans(self, inxscans,motionspec)
            if nargin==1,
                inxscans = 1:self.Nscan;
            end
            if max(inxscans)>self.Nscan
                error('Total scan: %d',self.Nscan);                
            end
                   
            for iscan = inxscans
                if nargin>2
                    MR_iscan = self.get_meanresp_insinglescan(iscan,motionspec);
                else
                    MR_iscan = self.get_meanresp_insinglescan(iscan);
                end
                MRfldnames = fieldnames(MR_iscan);     
                for ifldname =1:length(MRfldnames)
                    MRs(iscan).(MRfldnames{ifldname})=MR_iscan.(MRfldnames{ifldname});
                    %self.MRs(iscan).(MRfldnames{ifldname})=MR_iscan.(MRfldnames{ifldname});
                    %self.MRs(iscan)=setfield(self.MRs(iscan),MRfldnames{ifldname},MR_iscan.(MRfldnames{ifldname}));
                end                                
            end            
        end
        
        %% function [N C Evt CompleteCell inxsubset bspksucc] = collect_mresp_acrosscans(self, inxscans,outcell, bplot) 
        function [N C Evt CompleteCell inxsubset bspksucc] = collect_mresp_acrosscans(self, MRs,inxscans,outcell, bplot)            
            % function [N C CompleteCell bspksucc] = collect_mresp_acrosscans(self, inxscans, outcell, bplot)
            if nargin==1,
                inxscans = 1:self.Nscan;
            end            
            if max(inxscans)>self.Nscan
                error('Total scan: %d',self.Nscan);                
            end
            if nargin<3
                outcell=[];
            end
            
            dim_data = zeros(2,length(inxscans));
            
            inxi=1;            
            for iscan= inxscans
                mN = MRs(iscan).mN;
                [T nCell] = size(mN);                                
                dim_data(1,inxi) = T;
                dim_data(2,inxi) = nCell;
                inxi = inxi + 1;
            end
            
            % exclude cells with spike estimation failure
            min_nCell = min(dim_data(2,:));
            bspksucc = zeros(length(inxscans),min_nCell);
            inxi=1;            
            for iscan= inxscans
                mN = MRs(iscan).mN;
                bspk = mean(mN,1)>0;
                bspksucc(inxi,:) = bspk(1:min_nCell);
                inxi = inxi + 1;
            end
            cnt_spksucc = sum(bspksucc,1);
            
            %accumulated est spike data & calcium dF/F data
            CompleteCell =  find(cnt_spksucc==length(inxscans));
            CompleteCell = setdiff(CompleteCell, outcell);
            N = zeros(sum(dim_data(1,:)),length(CompleteCell));
            C = zeros(sum(dim_data(1,:)),length(CompleteCell));            
            Evt = zeros(sum(dim_data(1,:)),1);            
            inxsubset = [[0 cumsum(dim_data(1,1:end-1))]+1; cumsum(dim_data(1,:))];             
            inxi=1;              
            for iscan= inxscans
                mN = MRs(iscan).mN(:,CompleteCell);                    
                N(inxsubset(1,inxi):inxsubset(2,inxi),:)=mN;
                mC = MRs(iscan).mC(:,CompleteCell);                    
                C(inxsubset(1,inxi):inxsubset(2,inxi),:)=mC;  
                Evt(inxsubset(1,inxi):inxsubset(2,inxi))=MRs(iscan).events;  
                inxi = inxi + 1;
            end
            
            if nargin==4 && bplot,
                stepsize = 0.3;
                figure; hold on;
                chlabel=cell(1,size(N,2));
                for icell=1:size(N,2)
                    y=N(:,icell);
                    y=y(:)+(icell-1)*stepsize;    
                    plot(y(:))
                    chlabel{icell}= num2str(CompleteCell(icell));
                end
                for ii=1:length(inxscans)
                    plot(inxsubset(1,ii)*ones(1,size(N,2)+1),0:stepsize:stepsize*size(N,2),'r-.')
                end
                set(gca,'YTick',0:stepsize:(size(N,2)-1)*stepsize)
                
                set(gca,'YTickLabel',chlabel)
            end
                
        end
        
        %% function Nstim = get_frameresp_acrossscans(self, inxscans)
        function Nstim = get_frameresp_acrossscans(self, inxscans, motionspec)
            %  Nstim = get_frameresp_insinglescan(self,inxscan)
            % collect response in dFF and estimated spike in a single movie
            
         
            if max(inxscans)>self.Nscan
                error('Total scan: %d',self.Nscan);                
            end
            
            bquiet = true;
            for iscan = inxscans            
            
                if isfield(self.data(iscan).Params.stimparam,'stim1_samplesinNI')
                    stimtime = self.data(iscan).Params.stimparam.stim1_samplesinNI/self.data(iscan).Params.samplingfreq_NI;
                    blanktime = self.data(iscan).Params.stimparam.blank1_samplesinNI/self.data(iscan).Params.samplingfreq_NI;
                else
                    stimtime = self.data(iscan).Params.stimparam.stim_samplesinNI/self.data(iscan).Params.samplingfreq_NI;
                    blanktime = self.data(iscan).Params.stimparam.blank_samplesinNI/self.data(iscan).Params.samplingfreq_NI;
                end
                trialtime = stimtime+blanktime;
                nframe_trial = round(trialtime/self.data(iscan).Params.msperframe*1000);
                nCell=size(self.data(iscan).dFF,2);
                spec.frames = (-3:nframe_trial-1);
                spec.dataType = 'dFF';
                spec.nCell = nCell;
                if self.bawake
                    if nargin>2 
                        load(self.motionfn{iscan});
                        tmotionpix = sqrt(M(3,:).^2 + M(4,:).^2);
                        spec.motion.tmotion1 = tmotionpix;
                        spec.motion.motionthr1 = motionspec.motionthr1;
                        spec.motion.tmotion2 = [0 abs(diff(tmotionpix))];
                        spec.motion.motionthr2 = motionspec.motionthr2;
                        spec.motion.op='|';
                        spec.motion.bdisp = false;          
                    end
                end

                Xc = data_sort(self.data(iscan),spec,bquiet);
                clear spec;

                spec.frames=(-3:nframe_trial-1);
                spec.dataType='Nhat';
                spec.nCell=nCell;
                if self.bawake
                    if nargin>2 
                        load(self.motionfn{iscan});
                        tmotionpix = sqrt(M(3,:).^2 + M(4,:).^2);
                        spec.motion.tmotion1 = tmotionpix;
                        spec.motion.motionthr1 = motionspec.motionthr1;
                        spec.motion.tmotion2 = [0 abs(diff(tmotionpix))];
                        spec.motion.motionthr2 = motionspec.motionthr2;
                        spec.motion.op='|';
                        spec.motion.bdisp = false;          
                    end
                end
                [Xspk others]=data_sort(self.data(iscan),spec,bquiet);
                events = others.events{1};



                mXc1 = squeeze(mean(Xc{1},1));                
                smXc1 = std(mXc1,0,1);
                % cell selection procedure
                % select cell activity >0 in df/f and std(meanresp)< 3sigma
                learningvalid =zeros(1,nCell);
                for icell = 1:nCell
                    learningvalid(icell) = ~isnan(self.data(iscan).info.nhat_history(icell).posts(end));
                end
                inxincell = find(mean(mXc1,1)>0& smXc1<(mean(smXc1)+3*std(smXc1)) & learningvalid);
                inxoutcell = setdiff((1:nCell),inxincell);



                
                tf = (-3:nframe_trial-1)*self.data(iscan).Params.msperframe/1000;
                stimONframes = find(tf>=0 & tf<=stimtime);

                if isfield(self.data(iscan).Params.stimparam,'orientations')
                    ORI_deg = self.data(iscan).Params.stimparam.orientations;
                else
                    ORI_deg = self.data(iscan).Params.stimparam.orientation;                    
                end
                if isfield(self.data(iscan).Params.stimparam,'contrasts')
                    Cont = self.data(iscan).Params.stimparam.contrasts;
                else
                    Cont = self.data(iscan).Params.stimparam.contrast;
                end
                
                

                Xspk = Xspk{1}(:,stimONframes,:);                
                Xspk(:,:,inxoutcell)=0;

                trial_ori = self.data(iscan).Params.stimparam.cond.StimEventLUT(:,3);
                trial_cont = self.data(iscan).Params.stimparam.cond.StimEventLUT(:,1);


                Nstim(iscan).N_frame = Xspk;
                Nstim(iscan).ori = ORI_deg;
                Nstim(iscan).contrast = Cont;
%                 Nstim(iscan).trial_ori = trial_ori;
%                 Nstim(iscan).trial_cont = trial_cont;
                Nstim(iscan).timeTR = tf(stimONframes);
                Nstim(iscan).inxincell = inxincell;
                Nstim(iscan).events = events;
            end
            
        end %% end of get_frameresp_acrossscans
        
        %% function [N Evt CompleteCell inxsubset bspksucc] = collect_frameresp_acrosscans(self,Nstim, inxscans,outcell)      
        function [N Evt CompleteCell inxsubset bspksucc] = collect_frameresp_acrosscans(self,Nstim, inxscans,outcell)                        
            if nargin==1,
                inxscans = 1:self.Nscan;
            end            
            if max(inxscans)>self.Nscan
                error('Total scan: %d',self.Nscan);                
            end
            if nargin<3
                outcell={};
            end
            
            dim_data = zeros(2,length(inxscans));
            
            inxi=1;            
            for iscan= inxscans
                N_frame = Nstim(iscan).N_frame;
                [T frsize nCell] = size(N_frame);                                
                dim_data(1,inxi) = T;
                dim_data(2,inxi) = frsize;
                dim_data(3,inxi) = nCell;
                inxi = inxi + 1;
            end
            
            
            
            % exclude cells with spike estimation failure
            min_nCell = min(dim_data(3,:));
            bspksucc = zeros(length(inxscans),min_nCell);
            inxi=1;            
            for iscan= inxscans
                mN = squeeze(mean(Nstim(iscan).N_frame,2));
                bspk = mean(mN,1)>0;
                bspksucc(inxi,:) = bspk(1:min_nCell);
                inxi = inxi + 1;
            end
            cnt_spksucc = sum(bspksucc,1);
            
            %accumulated est spike data & calcium dF/F data
            CompleteCell =  find(cnt_spksucc==length(inxscans));
            CompleteCell = setdiff(CompleteCell, outcell);
            
                     
            N = zeros(dim_data(1,:)'*dim_data(2,:),length(CompleteCell));
            
            Evt = zeros(dim_data(1,:)'*dim_data(2,:),1);
            length_frameinscan = dim_data(1,:).*dim_data(2,:);
            
            inxsubset = [[0 cumsum(length_frameinscan(1:end-1))]+1; cumsum(length_frameinscan)];             
            inxi=1;              
            for iscan= inxscans
                Ni = Nstim(iscan).N_frame(:,:,CompleteCell); 
                Ni = permute(Ni,[2 1 3]);
                Ni = reshape(Ni, [size(Ni,1)*size(Ni,2), size(Ni,3)]);                
                N(inxsubset(1,inxi):inxsubset(2,inxi),:)=Ni;
                evt1 = repmat(Nstim(1).events,[1 dim_data(2,inxi)])';
                Evt(inxsubset(1,inxi):inxsubset(2,inxi)) = evt1(:);
                inxi = inxi + 1;
            end
            

                
        end
        
           
    end
end