function h = plotwithstimt(t,y, stimtimes,drawopts)
% function plotwithstimt(t,y, stimtimes,drawopts)
%
% t: time samples,
% y: signal(T x D); T= time course, D= # signal
% stimtimes: cell variable, each cell containing a stimulus onset and end time
% drawopts: other drawing options 

% Sangkyun Lee 2013-09-19

if length(t)~= size(y,1)
    error('# time sample must be the same as # y row');
end

for istim=1:length(stimtimes)
    
    stim=stimtimes{istim};
    minX = stim(1);
    maxX = stim(2);
    intv = (maxX-minX)/50;
    fX=[[minX:intv:maxX] [maxX:-intv:minX]];    
    
    minV = min(y(:));
    maxV = min(y(:))+ 1*(max(y(:))-min(y(:)));
    
    fY=[maxV*ones(1,length(fX)/2) minV*ones(1,length(fX)/2)];
    set(fill(fX,fY,[0.7 0.7 0.7]),'EdgeColor',[0.7 0.7 0.7]); hold on
end

if nargin<4,
    h=plot(t,y);
else
    cmdstring ='plot(t,y';
    for iopt=1:length(drawopts)
        if ~ isnumeric(drawopts{iopt})
            cmdstring = [cmdstring sprintf(',''%s''', drawopts{iopt})];
        else
            cmdstring = [cmdstring ', ' num2str(drawopts{iopt})];
        end
    end
    cmdstring = [cmdstring ')']
    h = eval(cmdstring); 
end