function s = setTimeStamp(s,timestamps)
%set the timestamps field in trial array s. The timestamps field will be acessed
%by sortStimEvent.
%
%make up the fields.
s.nevData = struct('TimeSpan',0,'TimeStampResolution',0,...
     'neurons',[],'events',[],'waves',[],'contvars',[]...
     );
s.nevFile = 'dummy';
s.nevFolder = 'c:\';

%assume timestamp is one dimension vector and set it into row format
if size(timestamps,1) > 1
    timestamps = timestamps';
end
    
s.nevData.events{1}.timestamps = timestamps;

