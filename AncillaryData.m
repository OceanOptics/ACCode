%% AncillaryData Abstract Interface Class
% 
% This object provides an abstract interface to multiple kinds of data
% types that might need to be used to process ACS data

classdef (Abstract) AncillaryData < handle
    
   properties (Abstract)
        Name
        Type
        Units
        % a timeseries object
        DataObject
        
        % preprocessing data
        SmoothData
% we dont' need this for tsg data
        SamplingFreq
        PPTimespan
        TransStartData;
        TransStartTime;
        TransEndData;
        TransEndTime;
               levelsMap;
   end
   properties

   end
   properties (Access = private)
       L   % logger
%         levelsMap;   % map to processing levels
       runningFSWmedian
       runningTSWmedian
        
    end
   methods (Abstract)

       qaqc( obj )
       bin( obj )
   end
   
   methods
       
       % try to define constructor here
       function obj = AncillaryData(nameIn, dataValuesIn, timestampsIn, unitsIn )
           
            obj.L = log4m.getLogger();
            obj.L.debug('AncillaryData.AncillaryData()','Created object');
           
            %disp('CALL ANCILLARY DATA CONSTRUCTOR')
            if nargin > 0 
                %disp('AncillaryData, args in')
                if ischar( nameIn )
                    obj.Name = nameIn;
                else
                    error('Name must be a string')
                end
                
                % check dataValues and timeStamps are valid first!
                obj.DataObject = timeseries(dataValuesIn, timestampsIn, 'name', nameIn);
                
                % assign units
                obj.Units = unitsIn;
                obj.DataObject.DataInfo.Unit = unitsIn;
            else
                % implicit call to superclass constructor?
                %disp('AncillaryData, no args in')
                obj.Name = '';
                obj.Type = '';
                
            end   % end if nargin > 0
            
            % define Quality codes 'lookup table' by calling defineQualityCode from
            % parent class
            %disp('set codes in AncillaryData')
            defineQualityCode(obj);
            
            % assign an initial quality code to each data record in matrix
            assignInitQualityCode(obj);
            
            % set up processing levels map
            keySet = {'raw', ...   %L1
                'preprocessed', ...%L2
                'binned', ...      %L3
                'filtered', ...    %L4
                'particulate', ... %L5
                'unsmoothed', ...  %L6
                'below750', ...    %L7
                'matchedWL', ...   %L8
                'corrected'};      %L9
            valueSet = {'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8', 'L9'};
            obj.levelsMap = containers.Map(keySet, valueSet);
            
        end   % end constructor
       
       function defineQualityCode(obj)
           %disp('call defineQualityCode() in AncillaryData')
           % Set the quality codes for AC Data Processing
           obj.DataObject.QualityInfo.Code = [ 1 2 3 4 9 ];
           obj.DataObject.QualityInfo.Description = {'Good', 'Unknown', 'Suspect', 'Bad', 'Missing'};
           
       end
              
       function assignInitQualityCode(obj)
           %disp('call assignInitQualityCode() in AncillaryData')
           % create a vector of ones ('Good') -- assume all data is good?
           obj.DataObject.Quality = ones( size(obj.DataObject.Data) );
       end
        
       function nRows = getSize(obj)
            nRows = obj.DataObject.Data.getdatasamplesize();
       end
        
       function obj = addData(obj, nameIn2, dataValuesIn2, timestampsIn2 )
           
            obj.L.debug('AncillaryData.addData', ...
                sprintf('timestampsIn 1: %s', datestr(timestampsIn2(1), 'HH:MM:SS:FFF')));
            obj.L.debug('AncillaryData.addData', ...
                sprintf('timestampsIn end: %s', datestr(timestampsIn2(end), 'HH:MM:SS:FFF')));
           obj.L.debug('AncillaryData.addData', ...
               sprintf('timestampsIn min: %s', datestr(min(timestampsIn2), 'HH:MM:SS:FFF')));
           obj.L.debug('AncillaryData.addData', ...
               sprintf('timestampsIn max: %s', datestr(max(timestampsIn2), 'HH:MM:SS:FFF')));
           
            % create new timeseries object from data coming in:
            % name is going to get wiped here as soon as we call append()
            tsIn = timeseries(dataValuesIn2, timestampsIn2, 'name', nameIn2);
            
            % assign initial quality codes
            tsIn.Quality = ones(size(dataValuesIn2));
            
           obj.L.debug('AncillaryData.addData', ...
               sprintf('Min Timestamp in new timeseries object: %s', datestr(min(tsIn.Time), 'HH:MM:SS:FFF')));
           obj.L.debug('AncillaryData.addData', ...
                sprintf('Max Timestamp in new timeseries object: %s', datestr(max(tsIn.Time), 'HH:MM:SS:FFF')));

            tsExist = obj.DataObject;
            
           obj.L.debug('AncillaryData.addData', ...
               sprintf('First timestamp in existing ts object: %s', datestr(tsExist.Time(1), 'HH:MM:SS:FFF')));
           
           obj.L.debug('AncillaryData.addData', ...
               sprintf('Last timestamp in existing ts object: %s', datestr(tsExist.Time(end), 'HH:MM:SS:FFF')));

     
            % append next time series to this one
            tsMerged = append(obj.DataObject, tsIn);
            
            % assign metadata from first ts to this new one. 
            % append() doesn't carry info across unless it matches exactly
            % and it doesn't carry .Name even if it DOES match.
            tsMerged.Name = obj.DataObject.Name;
            tsMerged.QualityInfo = obj.DataObject.QualityInfo;
            tsMerged.DataInfo = obj.DataObject.DataInfo;
     
            % reassign this DataObject to be the new ts3
            obj.DataObject = tsMerged;
            
            obj.L.debug('AncillaryData.addData', 'end of function');

       end
        
        function getInfo(obj)
            clear ts;
            ts = obj.DataObject;
            disp('Min Timestamp: getInfo')
            datestr(min(ts.Time))
            disp('Max Timestamp:getInfo')
            datestr(max(ts.Time))
            disp('Number of records:')
            length(ts.Time)
            disp('DataInfo')
            ts.DataInfo
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
        
        function plotData(obj)
            
            obj.L.debug('AncillaryData.plotData', 'Start');
            ts = obj.DataObject;
            scatter(obj.DataObject.Time, obj.DataObject.Data, '.')
            dynamicDateTicks();
            
        end   % end plotData
   end
     
end