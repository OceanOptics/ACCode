%% AncillaryInlineData Abstract Interface Class
% 
% This object provides an abstract interface to data collected "inline" 
% that needs to be processed in a similar way because of TSW/FSW
% transitions
% May 2016 - Added to code base; moved setting of levels map to here

classdef AncillaryInlineData < AncillaryData
   properties (Abstract)
        Name
        Type
        Units
        % a timeseries object
        DataObject
        
        % preprocessing data
        SmoothData
        levelsMap;
   end
   properties (Access = private)
       L   % logger
       runningFSWmedian
       runningTSWmedian
        
   end
    methods
           function obj = AncillaryInlineData(nameIn, dataValuesIn, timestampsIn, unitsIn)
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
            obj.L.debug('AncilllaryInlineData.ACData()','Created object');
            %setInfo(obj, nameIn);
                        

            
           end
           
        function setSmoothData(obj)
            
            obj.L.debug('AncillaryData.setSmoothData', 'start');            
            obj.SmoothData = filtfilt( ones(1,100)/100, 1, obj.DataObject.Data);
            
        end
        
        function plotSmoothData(obj)
            
            obj.L.debug('AncillaryData.plotSmoothData', 'start');
            
            scatter(obj.DataObject.Time, obj.SmoothData, '.');
            xlabel('Timestamp');
            ylabel(strcat(obj.Name, ' Smooth'));
            title(obj.Name);
            dynamicDateTicks;
            
        end  
           
           
   end
end