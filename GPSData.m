%GPSData - Data object to hold all processing information for latitude and longitude data
%GPSData inherits from AncillaryData
%
% Syntax:    obj = GPSData(nameIn, dataValuesIn, timestampsIn, unitsIn); 
% Inputs:
%    nameIn - name to identify object, i.e. 'gps'
%    dataValuesIn - actual data values
%    timestampsIn - actual timestamps for data
%    unitsIn - units of measurement for this data
%
%
% Example: 
%    gd = GPSData('LatAndLon', LatLon, fullDate, obj.GPSUnits); 
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

classdef GPSData < AncillaryData
    properties
        Name
        Type = 'gps';
        Units
        
        DataObject   % a timeseries object

        % processing variables for finding filtered
        SmoothData
        SamplingFreq
        PPTimespan
        TransStartData;
        TransStartTime;
        TransEndData;
        TransEndTime;
        levelsMap;
        
        var % Binned Data Variables 
        
    end
    properties (Access = private)
        L;   % Logger
    end;
    methods
        % constructor calls superclass constructor
        function obj = GPSData(nameIn, dataValuesIn, timestampsIn, unitsIn)
            %%% Pre-initialization %%%
            % Any code not using output argument (obj)
            
            %%% Object Initialization %%%
            % Call superclass constructor before accessing object
            % This statment cannot be conditionalized
                
            obj = obj@AncillaryData(nameIn, dataValuesIn, timestampsIn, unitsIn);
            
            %%% Post-initialization %%%
            % Any code, including access to the object
            %setInfo(obj, nameIn);
            obj.L = log4m.getLogger();
            obj.L.debug('GPSData.GPSData()','Created object');
            
        end
        
        % ----------------------------------------------------------------
        % qa/qc
        % ----------------------------------------------------------------         
        function qaqc(obj)
            disp('Stub for qaqc():');
            obj.Name
        end   % end QAQC()
        
        % ----------------------------------------------------------------
        % bin
        % ----------------------------------------------------------------        
        % flagThresholdIn -- number of minutes to check for gaps in data
        function bin(obj, timestampsIn, paramsIn)
            
            obj.L.debug('GPSData.bin()', 'Start method');

            level = obj.levelsMap('binned');
            flagThreshold = paramsIn.GPS_GAP_THRESHOLD;
            binSize = paramsIn.BIN_SIZE;
            methodIn = paramsIn.TSG_BIN_METHOD;
            latData = obj.DataObject.Data(:,1);
            lonData = obj.DataObject.Data(:,2);
            lTime = obj.DataObject.Time;
            
            % there are two possibilities for interpolating GPS data:
            % linear interpolation to ac bins, if the gaps between data points are large,
            % or by using standard binning method
            if strcmpi(methodIn, 'interpolate')
                
                obj.var.(level).BinnedLatData = interp1(lTime, latData, timestampsIn, 'linear');
                obj.var.(level).BinnedLonData = interp1(lTime, lonData, timestampsIn, 'linear');
                obj.var.(level).BinnedTimestamps = timestampsIn;
                obj.var.(level).LatBinFlags = ones(size(obj.var.(level).BinnedTimestamps));
                obj.var.(level).LonBinFlags = ones(size(obj.var.(level).BinnedTimestamps));
                
            elseif strcmpi(methodIn, 'bin')
                data = obj.DataObject.Data;
                time = obj.DataObject.Time;

                % limit to AC bin size
                startTimestamp = timestampsIn(1);
                endTimestamp = timestampsIn(end);
                
                % check for timestamps before first timestamps or after
                % last timestamp
                removeTimestamps = time < startTimestamp | time > endTimestamp;
                time(removeTimestamps) = [];
                data(removeTimestamps,:) = [];
          
                [binFlags, obj.var.(level).BinnedTimestamps,  ...
                    obj.var.(level).NumberBins, obj.var.(level).BinIndexNumbers, ...
                    binMean, binSampleSize, binSTD ] = bin(data, time, binSize,  timestampsIn );
                
                obj.var.(level).LatBinSampleSize = binSampleSize(:,1);
                obj.var.(level).LonBinSampleSize = binSampleSize(:,2);
                obj.var.(level).LatBinSTD = binSTD(:,1);
                obj.var.(level).LonBinSTD = binSTD(:,2);

            end;
            
            % check data is ok
            % ------------------------------------------------------------
            % Mark bins that have data with a predefined gap (i.e. >10 minutes)

              % do first for Lat Data

             [lat_data_changed, lat_flags_changed] = ...
                 checkGPSDataForGaps(data(:,1), time, binMean(:,1), obj.var.(level).BinnedTimestamps, binFlags(:,1));
             obj.var.(level).LatBinFlags = lat_flags_changed; 
             obj.var.(level).BinnedLatData = lat_data_changed;
             
             % do next for Lon
 
             [lon_data_changed, lon_flags_changed] = ...
                 checkGPSDataForGaps(data(:,2), time, binMean(:,2), obj.var.(level).BinnedTimestamps,  binFlags(:,2));
             obj.var.(level).LonBinFlags = lon_flags_changed;
             obj.var.(level).BinnedLonData = lon_data_changed;
            
            
          end   % end Bin()
        
        function plotData(obj)
            ts = obj.DataObject;
            scatter(ts.Data(:,2), ts.Data(:,1), [], ts.Time);
            xlabel('Longitude')
            ylabel('Latitude')
            title(obj.Name)
        end
        function plotBinnedData(obj, fignumIn)
            figure(fignumIn)
            grid on;
            hold on;
            scatter(obj.var.L3.BinnedLonData, obj.var.L3.BinnedLatData, ...
                [], obj.var.L3.BinnedTimestamps)
            xlabel('Longitude')
            ylabel('Latitude')
            title(obj.Name)
        end;
        
    end   % end methods
end   % end classDef