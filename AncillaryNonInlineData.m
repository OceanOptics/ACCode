%% AncillaryNonInlineData Abstract Interface Class
% 
% This object provides an abstract interface to data collected "inline" 
% that needs to be processed in a similar way because of TSW/FSW
% transitions

classdef AncillaryNonInlineData < AncillaryData
   properties (Abstract)
        Name
        Type
        Units
        % a timeseries object
        DataObject
        var % Binned Data Variables
        levelsMap;
   end
   properties (Access = private)
       L   % logger
   end
   methods (Abstract)
       bin( obj )
   end
   methods
           function obj = AncillaryNonInlineData(nameIn, dataValuesIn, timestampsIn, unitsIn)
            %%% Pre-initialization %%%
            % Any code not using output argument (obj)
            
            %%% Object Initialization %%%
            % Call superclass constructor before accessing object
            % This statment cannot be conditionalized
                
            obj = obj@AncillaryData(nameIn, dataValuesIn, timestampsIn, unitsIn);
            
            %%% Post-initialization %%%
            % Any code, including access to the object
            % adding variables to timeseries object.
            obj.L = log4m.getLogger();
            obj.L.debug('AncillaryNonInlineData.AncillaryNonInlineData()','Created object');
            %setInfo(obj, nameIn);
           end
   end
end