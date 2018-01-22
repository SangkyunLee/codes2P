function Params = loadParams(Params)
% function [Params] = loadParams(Params)
%
% Params.files.mainpath='D:\packages\onlineanalysis2P\test_data\TSeries-08292017-0850-000'
% samplingrate = 7.203%hz
% Params.Channel.inx_frtrg = 1;
% Params.Channel.inx_photo = 2;
% Params.Channel.nch =3;
% Params.files.stim_fn1='D:\packages\onlineanalysis2P\test_data\Retinotopy_1p5son_5x5.mat'
% Params.scan_mode = 'VA_spiral'
% Params.Nframes = Nframe;
% Params.msperframe = 1000/samplingrate;




Ch = Params.Channel;
nch = Ch.nch;

if isfield(Params,'Nframes')
    Nframe = Params.Nframes;
else
    Nframe =[];
end

%% load daq data
[DAQfn, Params] = get_DAQfn(Params);
[head,daq] = read_csv(DAQfn,nch);
% tmp = load('data_5khz.mat');
% daq = [tmp.time tmp.data];

headx = cell2mat(head(3:2:end));
t = daq(:,1)/1000; % insecond
daqx = daq(:,2:end);
Params.samplingfreq_NI = 1/diff(t(1:2));
Params.timeNI = t;    % in second



    
%% frame start time
y = daqx(:,headx==Ch.inx_frtrg);
fT = FrameT(t,y);
tinxf = fT.get_tstmp(Nframe);
if isempty(Nframe)
    Params.Nframes = length(tinxf);
end
Params.tinxf_NI = tinxf(:,1); %timestamp index in DAQ samples.

%% get photodiode time and visual stimulus parameters

pdt = PDT(t,daqx(:,headx==Ch.inx_photo));
Params = loadVStimParams(Params, pdt);


