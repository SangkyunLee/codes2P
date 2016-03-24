function hfig = plot2withstimt(data,idata,opts)
% function plot2withstimt(data,idata,opts)
%
% INPUT:
%     data: structure array
%     idata: index of data
%     opts.dattype;
%     opts.stepy;
%     opts.celllist;
%     opts.timeinterval=[startpoint endpoint]; in sample
%     opts.tunit ='frame' or 'msec' or 'sec'
%
% 2013-10-24, Sangkyun Lee    
    
dattype = opts.dattype;
stepy = opts.stepy;
cellist = opts.cellist;
N = length(cellist);
T = opts.timeinterval;
timeunit =opts.tunit;

if isempty(idata),
    idata=1;
end

msperframe = data(idata).Params.msperframe;
switch lower(timeunit)
    case 'frame'
        tinx = T(1):T(2);
        timex = [T(1):T(2)];        
        timescale = 1;
    case 'sec'
        tinx = round([T(1) T(2)]/msperframe*1000);
        tinx = tinx(1):tinx(end);
        tinx = tinx(find(tinx>0));       
        timex = tinx*msperframe/1000;
        timescale = msperframe/1000;
    case 'msec'
        tinx = round([T(1) T(2)]/msperframe);
        tinx = tinx(1):tinx(end);
        tinx = tinx(find(tinx>0));      
        timex = tinx*msperframe;
        timescale = msperframe;
    otherwise
        error('non pre-specified unit');

end
        
sig=data(idata).(dattype);



inx=find(data(idata).time.frame_start);
a=data(idata).time.stimtime1(inx);
inxstart=find(diff(a)>0);
inxend=find(diff(a)<0);
stimlen = length(inxstart)
stimtimes=cell(1,stimlen);
for ie=1:stimlen    
    stimtimes{ie}=[inxstart(ie) inxend(ie)]*timescale;
end



hfig =figure;
offsets_ch = [N-1:-1:0]*stepy;
offsets_chmat = repmat(offsets_ch,tinx(end)-tinx(1)+1,1);
sig = sig(tinx,cellist) + offsets_chmat;
plotwithstimt(timex,sig,stimtimes)
xlim([timex(1) timex(end)]);
ylim([offsets_ch(N) offsets_ch(1)]);
set(gca, 'YTick',fliplr(offsets_ch));
Yticklabel={};
for ich=1:size(sig,2)
    Yticklabel{ich} = num2str(size(sig,2)-ich+1);
end
 set(gca,'YTickLabel',Yticklabel)
 xlabel(timeunit);
