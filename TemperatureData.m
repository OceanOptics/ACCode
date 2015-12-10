%TemperatureData - Data object to hold all processing information for temperature data
%TemperatureData inherits from AncillaryData
%
% Syntax:    obj = TemperatureData(nameIn, dataValuesIn, timestampsIn, unitsIn); 
% Inputs:
%    nameIn - name to identify object, i.e. 'gps'
%    dataValuesIn - actual data values
%    timestampsIn - actual timestamps for data
%    unitsIn - units of measurement for this data
%
%
% Example: 
%    td = TemperatureData('Temperature', [TemperatureInstr, TemperatureBoat], fullDate, obj.TemperatureUnits);
%
% Other m-files required: AncillaryData
% Subfunctions: none
% MAT-files required: none
%
% See also: AncillaryData, TSGFileLoader

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% May 2015; Last revision: 13-08-15

%------------- BEGIN CODE --------------

classdef TemperatureData < AncillaryData
    properties
        Name
        Type = 'temperature';
        Units = 'degrees C';
        % a timeseries object
        DataObject
        
        SmoothData
        
        SamplingFreq
        PPTimespan
        TransStartData;
        TransStartTime;
        TransEndData;
        TransEndTime;
        
        var % Binned Data Variables
        
        levelsMap;
         
    end
    properties (Access = private)
        L;   % Logger
    end;
    methods
        % ----------------------------------------------------------------
        % constructor
        % ----------------------------------------------------------------
        
        function obj = TemperatureData(nameIn, dataValuesIn, timestampsIn, unitsIn)
            %%% Pre-initialization %%%
            % Any code not using output argument (obj)
            
            %%% Object Initialization %%%
            % Call superclass constructor before accessing object
            % This statment cannot be conditionalized
                
            obj = obj@AncillaryData(nameIn, dataValuesIn, timestampsIn, unitsIn);
            
            %%% Post-initialization %%%
            % Any code, including access to the object
            obj.L = log4m.getLogger();
            obj.L.debug('TemperatureData.TemperatureData()','Created object');

        end
        
        % ----------------------------------------------------------------
        % qa/qc
        % ----------------------------------------------------------------     
        function qaqc(obj)
            obj.L.debug('TemperatureData.qaqc()','Started method');
        end   % end QAQC()
        
        % ----------------------------------------------------------------
        % bin
        % ----------------------------------------------------------------
        function bin(obj, timestampsIn, paramsIn)
            
%             obj.L.debug('TemperatureData.bin()', 'Start method');

            level = obj.levelsMap('binned');
            flagThreshold = paramsIn.GPS_GAP_THRESHOLD;
            binSize = paramsIn.BIN_SIZE;
            methodIn = paramsIn.TSG_BIN_METHOD;
                        
            labTempData = obj.DataObject.Data(:,1);
            boatTempData = obj.DataObject.Data(:,2);
            tTime = obj.DataObject.Time;
            
            % there are two possibilities for interpoLabTemping Temperature data:
            % linear interpoLabTempion to ac bins, if the gaps between data points are large,
            % or by using standard binning method
            if strcmpi(methodIn, 'interpolate')
                
                obj.var.(level).BinnedLabTempData = interp1(tTime, labTempData, timestampsIn, 'linear');
                obj.var.(level).BinnedBoatTempData = interp1(tTime, boatTempData, timestampsIn, 'linear');
                obj.var.(level).BinnedTimestamps = timestampsIn;
                obj.var.(level).LabTempBinFlags = ones(size(obj.var.(level).BinnedTimestamps));
                obj.var.(level).BoatTempBinFlags = ones(size(obj.var.(level).BinnedTimestamps));
                
            elseif strcmpi(methodIn, 'bin')
                data = obj.DataObject.Data;
                time = obj.DataObject.Time;
                
                % data will now be equal to or smaller than ac timestamps
                
%                   disp('size of timestamps now')
%                   size(time)
%                   datestr(time(1))
%                   datestr(time(end))

                [binFlags, obj.var.(level).BinnedTimestamps,  ...
                    obj.var.(level).NumberBins, obj.var.(level).BinIndexNumbers, ...
                    binMean, binSampleSize, binSTD ] = ...
                    bin(data, time, binSize,  timestampsIn); %timestampsIn, binSize );
                
                LabTempBinFlags = binFlags(:,1);
                BoatTempBinFlags = binFlags(:,2);
                BinnedLabTempData = binMean(:,1);
                BinnedBoatTempData = binMean(:,2);
                obj.var.(level).LabTempBinSampleSize = binSampleSize(:,1);
                obj.var.(level).BoatTempBinSampleSize = binSampleSize(:,2);
                obj.var.(level).LabTempBinSTD = binSTD(:,1);
                obj.var.(level).BoatTempBinSTD = binSTD(:,2);

            end;
            
            % check data is ok
            % ------------------------------------------------------------
            % Mark bins that have data with a predefined gap (i.e. >10 minutes)
            
            % do first for Lab Temp Data

             [lab_temp_data_changed, lab_temp_flags_changed] = ...
                 checkGPSDataForGaps(data(:,1), time, BinnedLabTempData, obj.var.(level).BinnedTimestamps, LabTempBinFlags);
             obj.var.(level).LabTempBinFlags = lab_temp_flags_changed; 
             obj.var.(level).BinnedLabTempData = lab_temp_data_changed;
 
             [boat_temp_data_changed, boat_temp_flags_changed] = ...
                 checkGPSDataForGaps(data(:,2), time, BinnedBoatTempData, obj.var.(level).BinnedTimestamps, BoatTempBinFlags);
             obj.var.(level).BoatTempBinFlags = boat_temp_flags_changed;
             obj.var.(level).BinnedBoatTempData = boat_temp_data_changed;
 
 
          end   % end Bin()
        
        function plotData(obj)
  
            %disp('Stub for plotData() in TemperatureData:');
            ts = obj.DataObject;
            plot(ts.Time, ts.Data(:,1), 'k.')
            hold on
            scatter(ts.Time, ts.Data(:,2), 'b.')
            hold on
            dynamicDateTicks;
            xlabel('Timestamp')
            ylabel('Temperature - Lab & Intake')
            legend('Lab Temp', 'Intake Temp');
            title(obj.Name)

        end
        function plotBinnedData(obj, fignumIn)
            
            figure(fignumIn)
            grid on;
            hold on;
            plot(obj.var.L3.BinnedTimestamps, obj.var.L3.BinnedLabTempData, 'b');
            plot(obj.var.L3.BinnedTimestamps, obj.var.L3.BinnedBoatTempData, 'g');
            xlabel('Timestamp')
            ylabel('Temperature - Lab & Intake')
            legend('Lab/Instrument Temp', 'Intake Temp');
        end;
      
    end   % end methods
end   % end classDef