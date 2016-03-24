% classdef TPdata
%     properties
%         mainpath;% the root directory for data
%         stimulationtype; % continuous, flash, spontaneous
%         animal_status;
%         sessions;
%         
%         
%         
%         function field = get(self,fieldName)
%         switch fieldName
%             case {'mainpath','stimulationtype','refreshRate','resolution','gammaTable','debug','scrPixelPitch','displaySize'}
%                 field = e.(fieldName);
%             otherwise
%                 error('no suich field!')
%         end
%         
%         
%     end
% end
% 
