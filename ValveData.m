%ValveData - Data object to hold the all processing information for valve data
%ValveData inherits from AncillaryData
%
% Syntax:  valvedata = ValveData(nameIn, dataValuesIn, timestampsIn, unitsIn)
% Inputs:
%    nameIn - Description
%    dataValuesIn - Description
%    timestampsIn - Description
%    unitsIn - 
%
%
% Example: 
%     vd = ValveData('Valve1', ValveState, fullDate, obj.ValveUnits); 
%
% Other m-files required: AncillaryData
% Subfunctions: none
% MAT-files required: none
%
% See also: AncillaryData, FlowFileLoader

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% May 2015; Last revision: 13-08-15

%------------- BEGIN CODE --------------

classdef ValveData < AncillaryData
    properties
        Name
        Type = 'valve';
        Units
        % a timeseries object
        DataObject
        
        % preprocessing variables
        SmoothData
        SamplingFreq
        PPTimespan          
        TransStartData;
        TransStartTime;
        TransEndData;
        TransEndTime;
        levelsMap;
    end
    methods
      
        function obj = ValveData(nameIn, dataValuesIn, timestampsIn, unitsIn)
            
            %%% Pre-initialization %%%
            % Any code not using output argument (obj)
            
            %%% Object Initialization %%%
            % Call superclass constructor before accessing object
            % This statment cannot be conditionalized
                
            obj = obj@AncillaryData(nameIn, dataValuesIn, timestampsIn, unitsIn);
            
            %%% Post-initialization %%%
            % Any code, including access to the object

        end

        function findTransitions(obj)
            
            valveIndex = logical(obj.DataObject.Data);
    
            diffArray = valveIndex(2:end) - valveIndex(1:end-1);
            diffArray(end+1) = NaN;
            startsIndex = find(diffArray < 0);
            endsIndex = find(diffArray > 0);
            
            obj.TransStartTime = obj.DataObject.Time(startsIndex)
            obj.TransEndTime = obj.DataObject.Time(endsIndex)
            obj.TransStartData = obj.DataObject.Data(startsIndex);
            obj.TransEndData = obj.DataObject.Data(endsIndex);
        end        
          
        function qaqc(obj)
            disp('Stub for qaqc():');
            obj.Name
        end   % end QAQC()
        
        function bin(obj)
            disp('Stub bin():');
            obj.Name
        end   % end Bin()
        
        function plotTransitions(obj)
            hold on
            grid on
            plot(obj.DataObject.Time, obj.DataObject.Data)
            plot(obj.TransStartTime, obj.TransStartData, 'g*')
            plot(obj.TransEndTime, obj.TransEndData, 'r*')
            legend('valve state', 'valve opened/start FSW', 'valve close/end FSW')
            dynamicDateTicks
            hold off
            grid off            
        end

        function plotData(obj)
            
            ts = obj.DataObject;
            plot(obj.DataObject.Time, obj.DataObject.Data)
            xlabel('Timestamps')
            ylabel('Valve State')
            dynamicDateTicks();
            
        end   % end plotData
        
    end   % end methods
end   % end classDef