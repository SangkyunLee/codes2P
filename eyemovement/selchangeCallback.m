function selchangeCallback(h,~)
    k = get(h,'Parent');
    M = get(k, 'UserData'); 
    hR = get(get(M.ui.hgp1,'SelectedObject'),'String');
    M.DISP.mode =hR;
    set(k,'UserData',M);
end