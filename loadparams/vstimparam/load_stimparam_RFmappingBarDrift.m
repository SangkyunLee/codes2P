function stimparam = load_stimparam_RFmappingBarDrift(STIM,DAQ,constfldlist,trialfldlist)

paramconstants = STIM.constants;
paramtrial = STIM.trials(end);
ScreenRefreshrate = paramconstants.refreshRate;
stimparam.ScreenRefreshrate =  ScreenRefreshrate;
sampf = DAQ.samplingfreq_NI;


stim_samp = paramtrial.stimFrames*sampf/ScreenRefreshrate;
blank_samp = paramtrial.blankFrames*sampf/ScreenRefreshrate;
Barcenters = paramtrial.Barcenters;

nbarloc = length(Barcenters);
total_samp_singlesweep = (stim_samp + blank_samp)*nbarloc + paramtrial.postStimulusTime*sampf;
nsweep = length(paramtrial.BarDriftDRI);

stimparam.Barcenters = Barcenters;
stimparam.stim_samplesinNI = stim_samp;
stimparam.blank_samplesinNI = blank_samp;
stimparam.total_samp_singlesweep = total_samp_singlesweep;
stimparam.repetitions  = nsweep;


if ~exist('constfldlist','var') 
    constfldlist ={'displaySize','stimCenterDeg','stimRadiusDeg'};
end
for ifld = 1: length(constfldlist)
    stimparam.(constfldlist{ifld}) = paramconstants.(constfldlist{ifld});
end

if ~exist('trialfldlist','var') 
    trialfldlist ={'BarDriftDRI','BarWidth',...
        'StepSize','tempoFreq','numMotSteps',...
        'CheckSF','postStimulusTime'};
end
for ifld = 1: length(trialfldlist)
    stimparam.(trialfldlist{ifld}) = paramtrial.(trialfldlist{ifld});
end


evtcond = gen_evtinx(stimparam);
stimparam.evtcond = evtcond(:);
stimparam.expType = 'RFmappingBarDrift';
end

function evtcond = gen_evtinx(stimparam)
    Barcenters = stimparam.Barcenters;
    nbarloc = length(Barcenters);
    nrep = stimparam.repetitions;
    evtcond = repmat((1:nbarloc)', [1 nrep]);
    if stimparam.postStimulusTime>0,
        evtcond = [evtcond; zeros(1,nrep)];
    end

end