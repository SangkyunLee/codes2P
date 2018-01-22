function stimparam = load_stimparam_GratingMappingExperiment2PhotonbySang(STIM,DAQ)
paramconstants = STIM.constants;

ScreenRefreshrate = paramconstants.refreshRate;
stimparam.ScreenRefreshrate =  paramconstants.refreshRate;
sampf = DAQ.samplingfreq_NI;
stim_samp = paramconstants.stimFrames*sampf/ScreenRefreshrate;
blank_samp = paramconstants.blankFrames*sampf/ScreenRefreshrate;
total_samp = paramconstants.stimulusTime*sampf;
nrep = round(total_samp/(stim_samp + blank_samp));

stimparam.stim_samplesinNI = stim_samp;
stimparam.blank_samplesinNI = blank_samp;
stimparam.total_samplesinNI = total_samp;
stimparam.repetitions  = nrep;
stimparam.expType = 'GratingMappingExperiment2PhotonbySang';