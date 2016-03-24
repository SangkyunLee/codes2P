% Code for presenting different colorintensities
% modified by Sangkyun Lee, 5/8/2014

mode = 'Calib' % for measure gamma function of monitor
% mode = 'Test'; % After calculation of gammatable, check up whether luminance is linearlized.


colormode = 'gray';
% colormode = 'coloronly';
% colormode = 'all';


screens_ID= Screen('Screens');
screenNumber =1% max(screens_ID);

[xr,yr]=Screen(whichScreen,'WindowSize');
oldRes = SetResolution(screenNumber, 800, 600,72);
oldLevel = Screen('Preference', 'VisualDebugLevel',3); %turn off PsychToolbox Welcome Sign


[window, windowRect] = Screen('Openwindow', screenNumber, 0);
Screen('ColorRange', window, 1);
ifi = Screen('GetFlipInterval', window);

if strcmp(mode,'Calib')
    oldClut = LoadIdentityClut(window);
elseif strcmp(mode,'Test')
    origGammaTable = Screen('ReadNormalizedGammaTable',window);
    Screen('LoadNormalizedGammaTable',window,newGammaTable);
else
    error('Wrong mode');
end


maxPriorityLevel = MaxPriority(window);
Priority(maxPriorityLevel);

WaitSecs(2);

win_w = (windowRect(3) - windowRect(1));
win_h = (windowRect(4) - windowRect(2));

red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];
white = [1 1 1];

if strcmp(colormode,'gray')
    color_mat = [white]; 
elseif strcmp(colormode,'coloronly')
    color_mat = [red; green; blue; ]; 
else
    color_mat = [red; green; blue; white]; 
end

increments = 0:0.05:1;
luminance = cell(size(color_mat,1),length(increments));


% 
% imSize = 100;
% 
% viewingDist_cm = 62;
% foV_horizontal_cm = 52; 
% foV_vertical_cm = 32;
% 
% 
% viewingDist_cm = 50;
% foV_horizontal_cm = 44; 
% foV_vertical_cm = 28;
% 
% stimPos_inDeg_horizontal = 14;%10;
% stimPos_inRadians_horizontal = stimPos_inDeg_horizontal * pi/180;
% stimPos_inCM_horizontal = tan(stimPos_inRadians_horizontal) * viewingDist_cm;
% 
% conv2pixel_h = win_w/foV_horizontal_cm;
% stimPos_inPixel_horizontal = win_w/2 - round(stimPos_inCM_horizontal * conv2pixel_h);
% 
% 
% stimPos_inDeg_vertical = 5;
% stimPos_inRadians_vertical = stimPos_inDeg_vertical * pi/180;
% stimPos_inCM_vertical = tan(stimPos_inRadians_vertical) * viewingDist_cm;
% 
% conv2pixel_v = win_h/foV_vertical_cm;
% stimPos_inPixel_vertical = win_h/2 + round(stimPos_inCM_vertical * conv2pixel_v);
% 
% position_rect = CenterRectOnPoint([0 0 imSize imSize], stimPos_inPixel_horizontal, stimPos_inPixel_vertical);
% 
% 
% KbName('UnifyKeyNames');
% deviceNumberOperator = min(GetKeyboardIndices);


HideCursor;

for color_ind = 1:size(color_mat,1)
    
    color = color_mat(color_ind,:);
    
    for increment_ind = 1:length(increments)
      
        fill_color = increments(increment_ind)*color;
        
        %Screen('FillRect', window, fill_color, position_rect);
        Screen('FillRect', window, fill_color);
        %Screen('FrameRect', window, white, position_rect); 
        Screen('Flip', window);
        
%         while 1
%             
%             keyIsDown = KbCheck(deviceNumberOperator);
%             if keyIsDown
%                 keyIsDown = 0;
%                 break;
%             end
%         end    
        
        luminance{color_ind,increment_ind} = input('Luminance: ');
        WaitSecs(1);
    end
end
    

Screen('Flip', window);
Screen('CloseAll');
SetResolution(screenNumber, oldRes); %set monitor's previous resolution again
RestoreCluts

Screen('Preference', 'VisualDebugLevel', oldLevel); %so that PsychToolbox Initialzation Message will appear next time
clear Screen %turn off PsychDebugWindowConfiguration

ShowCursor;

save Luminance_projector.mat luminance
    