function [Out, J0]=sort_evtorder(evtlist,orders)
% function [Out J0]=sort_evtorder(evtlist,orders)
% evtlist=[cons' oris']
% orders={'descend','ascend'}
% 2016-03-17, Sangkyun Lee

J = zeros(size(evtlist));
J0 = 1 : size(evtlist,1);
for ievt = size(evtlist,2) : -1 :1
    [~, J(:,ievt)]=sort(evtlist(:,ievt),orders{ievt});    
end
for ievt = 1 : size(evtlist,2)
    J0 =J(J0,ievt);
end

Out=evtlist(J0,:);
        
            
    