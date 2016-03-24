function [ timestamps, codes ] = sortStimEvent(s,event,LUT)
%sort the stimulus events by tuning variables, orientation,spatial freq,etc
%Input: 
%s - data struct from matLoader 
%    s.matData
%    s.nexData
%    s.nevData
%note: s is 1-sized only. s.nevData.events is 1-sized only and timestamps are ISI filtered.
%      (s needs to be trimmed from matLoader,i.e,call trimNEVData before call this function ) 
%
%event - struct array
%    event.type - eventType to be sorted. e.g,{'Orientation','SpatialFreq'}
%    event.string - logic string on how to select the events by comparing the index.
%                   '>0'   select all values in the event variable. 
%                   '==1', select all values with condition index =1
%    event.operator - logic evaluation between adjecent event.string.
%                   - event(1).operator='&&' operates on event(1).string
%                   and event(2).string. 
%                   - last element will be ignored.
%change: pass Stim-Event-LUT as input parameter. this reduces significant
%time of getting it inside the script.
%
%Output :
% 
%

if length(s)> 1; error('array size of struct s has to be 1'); end;

m = length(s.nexData.markers);
for i = 1 : m
    if strcmp(s.nexData.markers{i}.name,'StimMarker')
        marker = s.nexData.markers{i};
        break;
    end
end

%selection event types
for i = 1 : length(event)
    event_type{i} = event(i).type;
end

%number of encoding(tuning) variables (including DIOValue)
nvar = length(marker.values);
var = struct;
for i = 1 : nvar
    var(i).type = marker.values{i}.name;
    var(i).string = '>0'; %no discrimation
    var(i).operator = '&'; %default is AND.
    if ~isempty(strmatch(var(i).type,event_type,'exact'))
        var(i).string = event(i).string;
        var(i).operator = event(i).operator;
    end    
end

%number of values for each variable -- they are generated as of same length .
nval = length(marker.values{i}.strings);

% LUT passed as input
if ~exist('LUT','var') || isempty(LUT)
    %transform the encoding variable array into a 2-d matrix
    LUT = zeros(nval,nvar); %stim event lookup table. row for values,column for variable
    %LUT = LUT'; %faster in manipulation of array ? 30/300, 10%
    %improvement.
    for i = 1 : nvar
        LUT(:,i) = str2double((marker.values{i}.strings)); %column vector containing encoding index
    end
    %use 16-bit format for sample index
    LUT = int16(round(LUT));
    %LUT = LUT';
end

%apply logic operator on each column of table
for i = 1 : nvar 
    S = eval(['LUT(:,i)' var(i).string]);
    if i == 1 ;  
        W = S;  
    else
        W = eval(['W' var(i-1).operator 'S']);
    end % W - accumulated logic array.
end

%retrieve timestamps from nevData --- already ISIFiltered.
timestamps = s.nevData.events{1}.timestamps;
%check the number of events are consistent, make sure the ISI filter worked fine

% if length(timestamps) ~= nval
%     error('Unequal counts of stimulus event: recorded=%d,presented=%d!',length(timestamps),nval);
% end

%logic selection of event-timestamps. column vector
if length(W)>length(timestamps)
    W(length(timestamps)+1:length(W)) = [];
end
    
timestamps = timestamps(W);
%encoding indices in the list for each variable (events x variable) array.
codes = LUT(find(W),:);





    







