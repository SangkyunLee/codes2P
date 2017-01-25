
function init_segment(videofn,frame)
% frame=[1:100];


%% pick up the threshold to find the pupil in the selected frame
self = load_video(videofn,frame);
% 
baseimg = single(self.data(:,:,1));
par.parname={'Mean'};
par.op1={'>p'};
par.pmaps = single(self.data(:,:,1));
par.fnsave='THRtemp.mat';
validate_thresholds(baseimg,par);

Isize = size(self.data(:,:,1));
M.data.img1 = mean(single(self.data),3);
figName = sprintf('Identify subregions for pupil');
M.fig = figure('Color', 'w', 'Name', figName, ...
			   'Units', 'pixels', 'Position', [100 100 1.2*Isize(2) 1.4*Isize(1)]);
M.ui.imgAxes = subplot('Position', [.1 .25 .83 .714]);                      

M.ui.export = uicontrol('Style','pushbutton','Units','normalized',...
    'position',[.4 .05 .15 .07],...
    'String','SET','Callback',@setCallback,'Enable','on');
M.ui.xy_pos = uicontrol('Style','text','Units','normalized',...
    'position',[.2 .05 .15 .07],...
    'String','','Enable','off');

set(M.fig,'WindowButtonDownFcn',{@wbd});

axes(M.ui.imgAxes); cla  
imagesc(M.data.img1);
colormap('gray');
colorbar;

set(M.fig,'UserData',M);
end
%---------------
% opt1.winsize=10;
% opt1.thr = 5;
% opt1.bdisp =false;
% 
% SEL.thr = opt1.thr ;
% SEL.op ='>';
% InitPar = gen_manual_pupilmask(self.data,1, SEL);
% out1 = track_rawdata(self.data,opt1,InitPar);

function setCallback(h,~)
    k=get(h,'Parent');
    M = get(k,'UserData');
    startpoint = M.ui.mouse.startpoint;
    endpoint = M.ui.mouse.stackpoint(:,1);
    pos_rec=sort([startpoint endpoint],2);
    assignin('base','pos_rec',pos_rec);
    close(M.fig);

end



% ---------------------------
function wbd(h,~)
    M = get(h,'UserData');    
    [m, n] = size(M.data.img1);
    % disp('down')


    pos = round(get(M.ui.imgAxes,'CurrentPoint'));
    pos = pos(1,1:2);
    
    if ~(pos(1)<1 || pos(1)>n || pos(2)<1 || pos(2)>m)
        M.ui.mouse.startpoint = pos';
        M.ui.mouse.stackpoint = [];
        % set the new values for the WindowButtonMotionFcn and
        % WindowButtonUpFcn
        set(h,'WindowButtonMotionFcn',{@wbm})
        set(h,'WindowButtonUpFcn',{@wbu})
        set(h,'UserData',M);
    end
end

% ---------------------------
function wbm(h,~)
    % executes while the mouse moves

    M = get(h,'UserData'); 
    [m, n] = size(M.data.img1);




    pos = round(get(M.ui.imgAxes,'CurrentPoint'));
    pos = round(pos(1,1:2));
    if ~(pos(1)<1 || pos(1)>n || pos(2)<1 || pos(2)>m)
        %opts.hAx = M.ui.imgAxes;        
        %hImg = overlayImage2(M.data.img1, M.data.roiimg, opts);
        
        axes(M.ui.imgAxes); cla  
        imagesc(M.data.img1);
        colormap('gray');
        colorbar;
        %---------------------------

        startpoint = M.ui.mouse.startpoint;
        endpoint = pos';
        Xpos = [startpoint(1) endpoint(1)];
        Ypos = [startpoint(2) endpoint(2)];

        lineX = [Xpos(1) Xpos(1)];
        lineY = [Ypos(1) Ypos(2)];
        line(lineX,lineY,'Color','y');
        lineX = [Xpos(1) Xpos(2)];
        lineY = [Ypos(2) Ypos(2)];
        line(lineX,lineY,'Color','y');
        lineX = [Xpos(2) Xpos(2)];
        lineY = [Ypos(1) Ypos(2)];
        line(lineX,lineY,'Color','y');
        lineX = [Xpos(1) Xpos(2)];
        lineY = [Ypos(1) Ypos(1)];
        line(lineX,lineY,'Color','y');

        M.ui.mouse.stackpoint = [pos' M.ui.mouse.stackpoint];
        set(h,'UserData',M);
    end
end


% ---------------------------
function wbu(h,~)
% executes when the mouse button is released

   

    M = get(h,'UserData');   
%     set(h,'WindowButtonMotionFcn','')
%     set(h,'WindowButtonUpFcn','') 
    
   
    [m, n] = size(M.data.img1);

    pos = round(get(M.ui.imgAxes,'CurrentPoint'));
    pos = round(pos(1,1:2));

    if ~(pos(1)<1 || pos(1)>n || pos(2)<1 || pos(2)>m)
        axes(M.ui.imgAxes); cla  
        imagesc(M.data.img1);
        colormap('gray');
        colorbar;


        startpoint = M.ui.mouse.startpoint;
        endpoint = pos';
        pos_rec=sort([startpoint endpoint],2);
        text_str = sprintf('[x=%d,y=%d - %d,%d]',...
            pos_rec(1,1), pos_rec(2,1),...
            pos_rec(1,2), pos_rec(2,2));
        set(M.ui.xy_pos,'String',text_str);        
        
        M.ui.mouse.stackpoint = [pos' M.ui.mouse.stackpoint];
        set(h,'UserData',M);

        
        set(h,'WindowButtonMotionFcn','')
        set(h,'WindowButtonUpFcn','') 
        

        

    end
end