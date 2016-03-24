function y = get_orifun(x,params)
% function y = get_orifun(x,params)
%
% defaults=struct('x0',0,'unit','deg','y0',0,'y1',1,'y2',1,'k1',1,'k2',1);
% y = y0 + y1*exp(k1.*(cos(x-x0)-1))+ y2*exp(k2.*(cos(x-x0+pi)-1));
%
% 2014-01-12, written by Sangkyun Lee

defaults=struct('x0',0,'unit','deg','y0',0,'y1',1,'y2',1,'k1',1,'k2',1);
% defaults=struct('x0',0,'unit','deg','model','von_mises','y0',0,'y1',1,'y2',1,'k1',1,'k2',1);
if nargin<2
    params = defaults;
else
    if isstruct(params)
        fnames = fieldnames(defaults);
        for i=1:length(fnames)
            if ~isfield(params,fnames{i})
                params.(fnames{i}) = defaults.(fnames{i});
            end
        end
    else
        error('params must be a struct');
    end    
end
x0 = params.x0;
y0 = params.y0;
y1 = params.y1;
y2 = params.y2;
k1 = params.k1;
k2 = params.k2;

if strcmp(params.unit,'deg')
    x = x/180*pi;
    x0 = x0/180*pi;
end

y = y0 + y1*exp(k1.*(cos(x-x0)-1))+ y2*exp(k2.*(cos(x-x0+pi)-1));
