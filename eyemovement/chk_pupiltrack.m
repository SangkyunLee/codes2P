function chk_pupiltrack(videofn,parfn, frames)
% function chk_pupiltrack(self,par)

% frames=5001:6000;



load(parfn);
par.out=fitout.out;
% do not transmit big data
par.out.sPIX = cell(1,length(fitout.out.sPIX));

par.fit = fitout.fitInfo;
par.searchInfo = searchInfo;

Xr=par.searchInfo.Xr;
Yr=par.searchInfo.Yr;
Nfr = par.searchInfo.srchdim(3);
if nargin>2 
    if max(frames)>Nfr || min(frames)<1
        error('Total # frame is %d',Nfr);
    end
    VR = load_video(videofn,frames,Xr,Yr);
    
else
    frames = 1:Nfr;
    VR = load_video(videofn,[],Xr,Yr);
end

%--------------------------

Isize = par.searchInfo.srchdim(1:2);

M.DISP.frameL = frames; % frame list
M.DISP.Cframe= frames(1); % current frame

PARs= zeros(3,Nfr);
for ii=1:Nfr
    if ~isempty(par.fit(ii).par)        
        PARs(:,ii)=par.fit(ii).par;
    end
end
M.DISP.PARs = PARs(:,M.DISP.frameL);

M.DISP.mode='NONE';
% --------------------------

figName = sprintf('Check pupil tracking');
M.fig = figure('Color', 'w', 'Name', figName, ...
			   'Units', 'pixels', 'Position', [100 100 2*Isize(2) 2.1*Isize(1)]);
M.ui.imgAxes = subplot('Position', [.05 .5238 .5 .476]);          

M.ui.CMxy = subplot('Position', [.05 .3 .9 .15]);  hold on;     
M.ui.Psize = subplot('Position', [.05 .05 .9 .15]); hold on;

dim = [.85 .44 .1 .01];
M.ui.ha1 = annotation(M.fig,'textbox',dim,'String','','FitBoxToText','on');
dim = [.85 .19 .1 .01];
M.ui.ha2 = annotation(M.fig,'textbox',dim,'String','','FitBoxToText','on');


M.ui.hgp1 = uibuttongroup('Units','normalized','Position',[.6 .8 .15 .15],'parent',M.fig,'SelectionChangeFcn',@selchangeCallback);
M.ui.rad1 = uicontrol('Style','Radio','Units','normalized','FontSize',11,...
    'String','NONE','position',[.1 .7 .8 .3],...
    'parent',M.ui.hgp1,'HandleVisibility','off');
M.ui.rad2 = uicontrol('Style','Radio','Units','normalized','FontSize',11,...
    'String','PIXEL','position',[.1 .4 .8 .3],...
    'parent',M.ui.hgp1,'HandleVisibility','off');
M.ui.rad3 = uicontrol('Style','Radio','Units','normalized','FontSize',11,...
    'String','FIT','position',[.1 .1 .8 .3],...
    'parent',M.ui.hgp1,'HandleVisibility','off');


M.ui.Fcursor = mrvSlider([.6 .65 .35 .1], 'Frame', ...
    'Range', [frames(1) frames(end)], 'IntFlag', 1, 'Value', frames(1), ...
    'FontSize', 11, 'Color', 'w','Callback','updatetracer');

Fstr = sprintf('[%d:%d]',frames(1),frames(end)); 
M.ui.Fedit = uicontrol('Style','edit','Units','normalized',...
        'position',[.6 .55 .2 .05],...
        'String',Fstr,'Enable','on','FontSize', 11);

M.ui.REEST = uicontrol('Style','pushbutton','Units','normalized',...
        'position',[.83 .55 .15 .05],...
        'String','Re-EST','Enable','on','FontSize', 11,'Callback',{@reest});
    
M.ui.EXPORT = uicontrol('Style','pushbutton','Units','normalized',...
        'position',[.78 .9 .15 .05],...
        'String','EXPORT','Enable','on','FontSize', 11,'Callback',{@export,parfn,fitout});    

M.data=VR.data;
cleardata(VR);
M.par =par;
set(M.fig,'UserData',M);
updatetracer(M);






end

%-----------------------
function export(h,~,parfn,fitout)
    k = get(h,'Parent');
    M = get(k, 'UserData'); 
    par=M.par;
    bj = ~cellfun(@isempty,par.out.sPIX);
    fitout.out.sPIX(bj) =par.out.sPIX(bj);
    fitout.out.CM = par.out.CM;
    fitout.fitInfo =par.fit;
    searchInfo = par.searchInfo;
    fitout.bmanual = true;
    save(parfn,'fitout','searchInfo','-v7.3');
    fitInfo=fitout.fitInfo;
    parfn1 = strrep(parfn,'fitinfo','smallfitinfo');
    save(parfn1,'fitInfo','searchInfo','-v7.3');
end


%-----------------
% 
function selchangeCallback(h,~)
    k = get(h,'Parent');
    M = get(k, 'UserData'); 
    hR = get(get(M.ui.hgp1,'SelectedObject'),'String');
    M.DISP.mode =hR;
    set(k,'UserData',M);
    updatetracer(M);
end


function reest(h,~)
    
    k = get(h,'Parent');
    M = get(k, 'UserData'); 

    eval(['frameL=' get(M.ui.Fedit,'String') ';']);
    frameL = frameL - M.DISP.frameL(1)+1;
    assert(frameL(1)<=frameL(end),'frame order should be incremental');

    set(h,'Enable','off');
    set(M.ui.Fedit,'Enable','off');
    
    reest_gui(M.data(:,:,frameL),frameL,M.fig);
    set(h,'Enable','on');
    set(M.ui.Fedit,'Enable','on');
end
