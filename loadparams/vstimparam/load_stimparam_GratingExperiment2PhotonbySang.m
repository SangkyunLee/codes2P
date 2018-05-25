function stimparam = load_stimparam_GratingExperiment2PhotonbySang(STIM,DAQ,constfldlist,trialfldlist)
paramconstants = STIM.constants;
paramtrial = STIM.trials(end);
ScreenRefreshrate = paramconstants.refreshRate;
stimparam.ScreenRefreshrate =  ScreenRefreshrate;
sampf = DAQ.samplingfreq_NI;
stim_samp = paramtrial.stimFrames*sampf/ScreenRefreshrate;
blank_samp = paramtrial.blankFrames*sampf/ScreenRefreshrate;
total_samp = paramtrial.stimulusTime*sampf;
ntrial = round(total_samp/(stim_samp + blank_samp));

stimparam.stim_samplesinNI = stim_samp;
stimparam.blank_samplesinNI = blank_samp;
stimparam.total_samplesinNI = total_samp;
stimparam.repetitions  = ntrial;


if ~exist('constfldlist','var') 
    constfldlist ={'displaySize','stimCenterDeg','stimRadiusDeg',...
        'stimMarginDeg','gratingWaveform'};
end
for ifld = 1: length(constfldlist)
    stimparam.(constfldlist{ifld}) = paramconstants.(constfldlist{ifld});
end

if ~exist('trialfldlist','var') 
    trialfldlist ={'contrast','orientation','spatialFreq','tempoFreq','direction','DIOValue'};
    stimparam.DIOcondseq =trialfldlist;
end
for ifld = 1: length(trialfldlist)
    stimparam.(trialfldlist{ifld}) = paramtrial.(trialfldlist{ifld});
end


stimparam.expType = 'GratingExperiment2PhotonbySang';