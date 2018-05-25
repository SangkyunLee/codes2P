function stimparam = load_stimparam_SquareMappingExperiment2Photon(STIM,DAQ,fldlist)
% stimparam = load_stimparam_SquareMappingExperiment2Photon(STIM,DAQ,fldlist)

    paramconstants = STIM.constants;
    paramtrial = STIM.trials(end);
    ScreenRefreshrate = paramconstants.refreshRate;
    stimparam.ScreenRefreshrate =  ScreenRefreshrate;
    sampf = DAQ.samplingfreq_NI;
    stim_samp = paramconstants.stimFrames*sampf/ScreenRefreshrate;
    blank_samp = paramconstants.blankFrames*sampf/ScreenRefreshrate;
    total_samp = paramconstants.stimulusTime*sampf;
    ntrial = round(total_samp/(stim_samp + blank_samp));

    stimparam.stim_samplesinNI = stim_samp;
    stimparam.blank_samplesinNI = blank_samp;
    stimparam.total_samplesinNI = total_samp;
    stimparam.repetitions  = ntrial;



    if ~exist('fldlist','var') 
        fldlist ={'validTrial','dotNumX','dotNumY','stimCenterX','stimCenterY','dotSize','rndfn','dotLocations','dotColors','direction','DIOValue'};
    end
    for ifld = 1: length(fldlist)
        if isfield(paramconstants,fldlist{ifld})
            stimparam.(fldlist{ifld}) = paramconstants.(fldlist{ifld});
        elseif isfield(paramtrial,fldlist{ifld})
            stimparam.(fldlist{ifld}) = paramtrial.(fldlist{ifld});
        end
    end
    LUT = condLUT(stimparam);
    stimparam.evtcond = LUT;
    stimparam.expType = 'SquareMappingExperiment2Photon';
end

function LUT = condLUT(VStimparam)

    n = cellfun(@size,VStimparam.dotLocations,'UniformOutput',false);
    n = cell2mat(n');
    nmaxdot = max(n(:,2));

    if nmaxdot==1,
        dotLoc = cell2mat(VStimparam.dotLocations);
        [xgrid, ~, xinx] = unique(dotLoc(1,:));
        [ygrid, ~, yinx] = unique(dotLoc(2,:));
        dotClr = cell2mat(VStimparam.dotColors);
        [dotclr, ~, cinx] = unique(dotClr);
        nx = length(xgrid);
        ny = length(ygrid);
        nc = length(dotclr);

        LUT =[xinx yinx cinx sub2ind([nx ny nc],xinx,yinx,cinx)];  
    else
        error('not implemented yet for simultaneouse multi-dot stimulation');
    end
end