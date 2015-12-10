%SalinityData - Data object to hold all processing information for salinity data
%SalinityData inherits from AncillaryData
%
% Syntax:   obj = SalinityData(nameIn, dataValuesIn, timestampsIn, unitsIn)
% Inputs:
%    nameIn - name to identify object, i.e. 'salinity'
%    dataValuesIn - actual data values
%    timestampsIn - actual timestamps for data
%    unitsIn - units of measurement for this data
%
%
% Example: 
%    sd = SalinityData('Salinity', Salinity, fullDate, obj.SalinityUnits); 
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

classdef SalinityData < AncillaryData
    properties
        Name            % name for printing on plots, etc.
        Type = 'gps';   % type of file this data comes from
        Units           % units of measurement
        
        DataObject      % a timeseries object that will hold data, time
        
        SmoothData      
        SamplingFreq
        PPTimespan
        TransStartData;
        TransStartTime;
        TransEndData;
        TransEndTime;
        levelsMap;
        
        
        var            % Binned Data Variables 
    end
    properties (Access = private)
        L;   % Logger
    end;
    methods
        
        % constructor
        function obj = SalinityData(nameIn, dataValuesIn, timestampsIn, unitsIn)
            
            %%% Pre-initialization %%%
            % Any code not using output argument (obj)
            
            %%% Object Initialization %%%
            % Call superclass constructor before accessing object
            % This statment cannot be conditionalized
                
            obj = obj@AncillaryData(nameIn, dataValuesIn, timestampsIn, unitsIn);
            
            %%% Post-initialization %%%
            % Any code, including access to the object
            obj.L = log4m.getLogger();
            obj.L.debug('SalinityData.SalinityData()','Created object');
        end

        function qaqc(obj)
            obj.L.debug('SalinityData.qaqc()', 'Start method');
        end   % end QAQC()
        
        % ----------------------------------------------------------------
        % bin
        % ----------------------------------------------------------------
        % flagThresholdIn -- number of minutes to check for gaps in data
        function bin(obj, timestampsIn, paramsIn)
            
            obj.L.debug('SalinityData.bin()', 'Start method');

            level = obj.levelsMap('binned');
            flagThreshold = paramsIn.GPS_GAP_THRESHOLD;
            binSize = paramsIn.BIN_SIZE;
            methodIn = paramsIn.TSG_BIN_METHOD;
            
            salData = obj.DataObject.Data(:,1);
            salTime = obj.DataObject.Time;
            
            % there are two possibilities for interpolating GPS data:
            % linear interpolation to ac bins, if the gaps between data points are large,
            % or by using standard binning method
            if strcmpi(methodIn, 'interpolate')
                
                obj.var.(level).BinnedData = interp1(salTime, salData, timestampsIn, 'linear');
                obj.var.(level).BinnedTimestamps = timestampsIn;
                obj.var.(level).BinFlags = ones(size(obj.var.(level).BinnedTimestamps));
                
            elseif strcmpi(methodIn, 'bin')
                data = obj.DataObject.Data;
                time = obj.DataObject.Time;
                size(time)
                size(data)
                % limit to AC bin size
                startTimestamp = timestampsIn(1);
                endTimestamp = timestampsIn(end);
                
                datestr(startTimestamp)
                datestr(endTimestamp)
                
                % check for timestamps before first timestamps or after
                % last timestamp
                removeTimestamps = time < startTimestamp | time > endTimestamp;
                time(removeTimestamps) = [];
                data(removeTimestamps,:) = [];
                
%                 size(time)
%                 size(data)
          
                [binFlags, obj.var.(level).BinnedTimestamps,  ...
                    obj.var.(level).NumberBins, obj.var.(level).BinIndexNumbers, ...
                    binMean, binSampleSize, binSTD ] = bin(data, time, binSize,  timestampsIn );
                
                obj.var.(level).BinSampleSize = binSampleSize(:,1);
                obj.var.(level).BinSTD = binSTD(:,1);

            end;
            
             [data_checked, flags_checked] = ...
                 checkGPSDataForGaps(data, time,  binMean(:,1), obj.var.(level).BinnedTimestamps, binFlags(:,1));
             obj.var.(level).BinFlags = flags_checked; 
             obj.var.(level).BinnedData = data_checked;
          
          end   % end Bin()
          
          
          function plotData(obj)
  
            sal = obj.DataObject;
            scatter(sal.Time, sal.Data, 'b.')
            hold on
            dynamicDateTicks;
            xlabel('Timestamp')
            ylabel('Salinity')
            title(obj.Name)

        end
          
        % ---------------------------------------------------------------
        % plotBinnedData()
        % ---------------------------------------------------------------
        function plotBinnedData(obj, fignumIn)
            figure(fignumIn)
            grid on;
            hold on;
            plot(obj.var.L3.BinnedTimestamps, obj.var.L3.BinnedData, 'b');
            xlabel('Binned Timestamps');
            ylabel('Binned Data')
            title('Binned Salinity Data')
        end;
    end   % end methods
end   % end classDef