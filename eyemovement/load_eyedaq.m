function [daqdata,vidtime,Info] = load_eyedaq(eye_dir,recordnum)
% function [daqdata,vidtime,Info] = load_eyedaq(eye_dir,recordnum)
% vidtime is the relative time from the daqtriggered time
% 2016-04-05 Sangkyun Lee

    dlogfn = fullfile(eye_dir,sprintf('%04d.dlog',recordnum));
    vlogfn = fullfile(eye_dir,sprintf('%04d.vlog',recordnum));
    load(fullfile(eye_dir,sprintf('%04d.mat',recordnum)))


    vfid =fopen(vlogfn,'rt');
    vlog=single(fscanf(vfid,'%d, %f, %f\n'));
    fclose(vfid);
    vlog = reshape(vlog,[3 length(vlog)/3])';


    dfid = fopen(dlogfn,'rb');
    dlog = fread(dfid,'double');
    fclose(dfid);
    DAQTriggerTime=dlog(2);

    dlog =dlog(4:end);
    Nch=2;
    daqdata=reshape(dlog,[Nch length(dlog)/Nch])';

    Info.daqtriggertime = datestr(DAQTriggerTime,'mmmm dd, yyyy HH:MM:SS.FFF AM');
    Info.videotriggertime = datestr(datenum(data_params.video.InitialTriggerTime),'mmmm dd, yyyy HH:MM:SS.FFF AM');

    % datestr(TriggerTime, 'HH:MM:SS.FFF')
    % datestr(datenum(data_params.video.InitialTriggerTime), 'HH:MM:SS.FFF')
    trigger_offset = str2double(datestr(DAQTriggerTime - datenum(data_params.video.InitialTriggerTime),'SS.FFF'));
    %daqtime = daqdata(:,1);
    vidtime = vlog(:,2)-trigger_offset;

end
% %------ find the timestamp for the first visual stimulation
% 
% daqfr = data_params.daqinfo.Rate;
% sig = daqdata(:,2); 







