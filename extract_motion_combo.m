function [speed, dist_twin, IMG_motion, outparam] = extract_motion_combo(DAQpath,FMOTpath,pixel_resolution, frame_period, twin,channel, bfig)

% rotation pattern (clockwise or counter clockwise)
ROT1 = [0 0 0 1 1 1 1 0];
ROT2 = [1 0 1 1 0 1 0 0];
ROTAROADSIZE_RAD = 8; % in cm
ENCODER_RESOLUTION = 8000; % MAXIMUM CYCLE 


DAQfstr = dir([DAQpath filesep '*.csv'])
if isstruct(DAQfstr)
    DAQfn = DAQfstr.name;
else
    error('no csv file in the specified directory');
end

fullfn_DAQ = fullfile(DAQpath,DAQfn)
fid = fopen(fullfn_DAQ);
C1 = textscan(fid, '%s %s%d, %s%d, %s%d, %s%d\n');
DAQ_data = textscan(fid, '%f, %f, %f, %f, %f','CollectOutput', 1);
fclose(fid)    
DAQ_data = DAQ_data{1};    
DAQ_data(:,1) = DAQ_data(:,1)/1000;
DAQ_sampleperiod =  1/diff(DAQ_data(1:2,1));
% figure; plot(DAQ_data(:,1),DAQ_data(:,4))


Ymov=[]; Xmov=[];
diffPhase=[];
flist = dir([FMOTpath filesep 'motionparameter*.mat']);
if ~isempty(flist)
    for iseg=1:length(flist)
        motionfn=flist(iseg).name;
        fullfn_motion = fullfile(FMOTpath,motionfn);
        a=load(fullfn_motion)
        diffPhase=[diffPhase a.M(2,:)];
        Ymov=[Ymov a.M(3,:)];
        Xmov=[Xmov a.M(4,:)];
    end
else
    error('no motionparameter*.mat files in the specified directory');
end
if length(find(isnan(diffPhase)))>1
    error('motion parameters include nan');
end


len1 = length(Ymov);


%%
ROT_motion = DAQ_data(:,[1 channel]);
% figure; plot(ROT_motion(:,1), ROT_motion(:,2:3));

IMG_motion(:,1) = [1:len1]*frame_period;
IMG_motion(:,2) = sqrt(Xmov.^2+Ymov.^2)*pixel_resolution;


A=ROT_motion(:,2)-mean(ROT_motion(:,2));
A=(sign(A/max(abs(A)))+1)/2;

B=ROT_motion(:,3)-mean(ROT_motion(:,3));
B=(sign(B/max(abs(B)))+1)/2;
T_rot = ROT_motion(:,1);

%  figure; plot(T_rot,[A B])

statechange_inxA=find(diff(A)~=0);
statechange_inxB=find(diff(B)~=0);
statechange_inx = sort([statechange_inxA;statechange_inxB]);

MA=[A(statechange_inxA) B(statechange_inxA) A(statechange_inxA+1)  B(statechange_inxA+1)];
MB=[A(statechange_inxB) B(statechange_inxB) A(statechange_inxB+1)  B(statechange_inxB+1)];

Time_statechange = T_rot(statechange_inx);
M=[A(statechange_inx) B(statechange_inx) ];



Maug = [M(1:end-3,:) M(2:end-2,:) M(3:end-1,:) M(4:end,:) ];
direction=[];
RotTime =Inf*ones(size(Maug,1),1);
for istate = 1: size(Maug,1)
    Mi = Maug(istate,:);
    for ishift = 0 : 2 : 6
        ROT1sh= circshift(ROT1',ishift)';
        ROT2sh= circshift(ROT2',ishift)';
        if sum(abs(Mi-ROT1sh))==0,
            direction(istate)=1;
            RotTime(istate)=Time_statechange(istate+3) - Time_statechange(istate);
        end
        if sum(abs(Mi-ROT2sh))==0,
            direction(istate)=-1;
            RotTime(istate)=Time_statechange(istate+3) - Time_statechange(istate);
        end
    end
end


speed=(ROTAROADSIZE_RAD*2*pi)./(RotTime*ENCODER_RESOLUTION);
len=min(length(speed),length(direction))
%  figure; 
 
% close all
if exist('bfig','var') && bfig
    figure; hold on;
    plot(ROT_motion(:,1),ROT_motion(:,2:3))
    plot(IMG_motion(:,1),IMG_motion(:,2),'k')
    plot(Time_statechange(2:len+1), speed(1:len).*direction(1:len)','y','LineWidth',2)
    plot(Time_statechange(2:len+1), speed(1:len),'r','LineWidth',1.5)
end


Time_statechange = Time_statechange(2:len+1);
speed_tstamp = speed(1:len);
direction = direction(1:len);
statechange_inx = statechange_inx(2:len+1);
speed = zeros(size(DAQ_data,1),1);
speed(statechange_inx)=speed_tstamp;

outparam.Time_statechange=Time_statechange;
outparam.speed_tstamp = speed_tstamp;
outparam.direction = direction;
outparam.statechange_inx = statechange_inx;
IMG_motion = IMG_motion(:,2);
 
if exist('twin','var')
    %calculation the distance in cm during the given time twin
    Ntwin = round(DAQ_sampleperiod*twin);
    speed_twin = ones(Ntwin,1)*twin;
    dist_twin = conv(speed,speed_twin);
    dist_twin = dist_twin(round(Ntwin/2)+1:length(speed)+round(Ntwin/2));

    if exist('bfig','var')
        plot(DAQ_data(:,1),dist_twin,'c','LineWidth',2)
        ylim([0 10])
    end
 else
     dist_twin=[];
 end
 
 
 
 return;
 

 
 
             
           

