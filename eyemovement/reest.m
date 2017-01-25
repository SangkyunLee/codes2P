% function reest(h,~,VR)
%     
%     k = get(h,'Parent');
%     M = get(k, 'UserData'); 
%     eval(['frameL=' get(M.ui.Fedit,'String') ';']);
%     set(h,'Enable','off');
%     set(M.ui.Fedit,'Enable','off');
%     
%     frameL = frameL - M.DISP.frameL(1)+1;
%     reest_gui(VR,frameL,M);
%     set(h,'Enable','on');
%     set(M.ui.Fedit,'Enable','on');
% end