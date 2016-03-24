function out = gen_manual_pupilmask(MOV, frame, SEL)
% function out = gen_manual_pupilmask(MOV, frame, SEL)
% SEL.thr = 30;
% SEL.op = '>'; ==> thr>Img
% 2016-03-22, Sangkyun Lee

assert(isstruct(SEL),'SEL should be a struct with fields ''thr'' and ''op''');
switch SEL.op
    case {'>','<','<=','>='}
    otherwise
        error('not proper operator');
end

fImg = single(MOV(:,:,frame));
fmap = zeros(size(fImg));
thr = SEL.thr;
cmd = sprintf('fmap(thr%sfImg)=1;',SEL.op);
eval(cmd);
L1 = bwlabel(double(fmap),8);
hfig1 = figure; imagesc(L1);
title('SELECT PUPIL');
uiwait(hfig1);
prompt = {'Enter pupil cluster ID'};
dlg_title = 'PUPIL ID';
num_lines = 1;
def = {'1'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
ID = str2double(answer{1});
L1(L1(:)~=ID)=0;
Cpix = find(L1(:)==ID);
[X0, Y0]=cmass(L1);

out.CM=[X0,Y0];
out.pix = Cpix;
