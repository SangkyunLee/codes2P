function [stimstart, saminx] = identify_stimstart(daqtime, PDsig)
%function [stimstart, saminx] = identify_stimstart(daqtime, PDsig)
% stimstart : time in daqtime unit
% sameple inx 
% 2016-04-05 Sangkyun Lee

len=round(length(PDsig)*0.2);
h=figure; plot(daqtime(1:len,1),PDsig(1:len)); xlabel('time (sec)')
waitfor(h)


prompt = {'Threshold','Time start'};
dlg_title = 'Input a starting point for the photo-diode';
num_lines = 2;
def = {'2','5'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
thresh_PD = eval(answer{1});
inxt = eval(answer{2});

saminx =inxt;

while PDsig(saminx)<thresh_PD
    saminx=saminx+1;
    if mod(saminx,100)==1, pause(0.01); end
end


stimstart = saminx* diff(daqtime(1:2));
end