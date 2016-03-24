function [beta r handles] = plot_linearreg(X,hax, plotopt)
% function [beta r] = plot_linearreg(X,hax)    
    
    if nargin<2 || isempty(hax)
        hs = plot(X(:,1),X(:,2),'.');
    else
       hs = plot(hax,X(:,1),X(:,2),'.');
    end
    
    if nargin==3  && isfield(plotopt,'scatter'),
        
        scatteropt = plotopt.scatter;
        Sfields = fieldnames(scatteropt);
        for iprop = 1 : length(Sfields)
            set(hs,Sfields{iprop},scatteropt.(Sfields{iprop}));
        end
    end
    handles.hs=hs;
    
        
        
    Xv=[X(:,1) ones(length(X(:,1)),1)];
    Y=X(:,2);
    beta=Xv\Y;
    r=corr([Xv(:,1), Y]);
    xline=[min(X(:,1)):(max(X(:,1))-min(X(:,1)))/10:max(X(:,1))];
    hold on; 
    if nargin<3  || (~isfield(plotopt,'line') || ~isfield(plotopt.line,'Color'))
        hl=plot(xline,beta(1)*xline+beta(2),'r');
    else
        hl=plot(xline,beta(1)*xline+beta(2));
    end
    str=sprintf('Y=%.2fX+%.2f; R=%.2f',beta(1),beta(2),r(1,2));
    if nargin==3  && isfield(plotopt,'line'),
        lineopt = plotopt.line;
        Lfields = fieldnames(lineopt);
        for iprop = 1 : length(Lfields)
            set(hl,Lfields{iprop},lineopt.(Lfields{iprop}));
        end        
    end
    handles.hl=hl;
    title(str);
end