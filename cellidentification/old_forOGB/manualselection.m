function [inx_inroiwithboundary,inxboundary]=manualselection(cellids, Img, inx_inroiwithboundary,inxboundary)
% Img=template_org;
for inxcell=cellids
    clear global clststruct
    global clststruct
    hFig=figure('Position',[100 678 560 420]);
    mansep(hFig,Img,inx_inroiwithboundary{inxcell},inxboundary{inxcell},inxcell);
    uiwait(hFig);%uiwait(gcf)
%     close(hFig);
    pause(0.1);
    % figure; imagesc(template_org)
    % a=template_org;%zeros(size(template_org));
    % selinx = clststruct.sel_boundvoxinFOV{1};%union(clststruct.sel_invoxinFOV{2},clststruct.sel_boundvoxinFOV{2});
    % a(selinx)=1000;
    % figure; imagesc(a)
    inx_inroiwithboundary{inxcell} = union(clststruct.sel_invoxinFOV{1},clststruct.sel_boundvoxinFOV{1})';
    inxboundary{inxcell} = clststruct.sel_boundvoxinFOV{1};
    inx_inroiwithboundary{length(inx_inroiwithboundary)+1} = union(clststruct.sel_invoxinFOV{2},clststruct.sel_boundvoxinFOV{2})';
    inxboundary{length(inxboundary)+1} = clststruct.sel_boundvoxinFOV{2};
    
end
inxempty=[];
temp1={};
temp2={};
count=1;
for inxcell=1:length(inxboundary)
    if ~isempty(inxboundary{inxcell})        
        a=inxboundary{inxcell};        
        temp1{count}=a(:)';
        b=inx_inroiwithboundary{inxcell};
        temp2{count}=b(:)';
        count = count+1;
    end
end
inxboundary=temp1;
inx_inroiwithboundary=temp2;