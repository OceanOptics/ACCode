% PROCESSEDDATA is a class that contains the data for AC Data
% Processing.  
%
%
% Other m-files required: spectral_unsmooth, ResidTempScatCorr,
% AttTempCorr, getPsiT
% Other files required:   'Sullivan_etal_2006_instrumentspecific.xls'
% Subfunctions: spectral_unsmooth, ResidTempScatCorr, AttTempCorr, getPsiT
% MAT-files required: none
%
% See also: spectral_unsmooth, ResidTempScatCorr, AttTempCorr, getPsiT

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% May 2015; Last revision: 13-08-15
% July 2016 - Major revision

% ------------- BEGIN CODE --------------

classdef ProcessedData < handle
    
    properties (Access = private)
        numWavelengths;   %# number of wavelengths for this ac meter
    end   % end private properties
    
    properties
        var               %# all processing variables
        meta              %# processing metadata
        levelsMap         %# a map of the level names and numbers
        levelsFlags       %# a set of flags showing which levels are in use
        cExists           %# boolean
        aExists           %# boolean
        L                 %# logger
    end   % end public properties
    
    methods
        
    %% ----------------------------------------------------------------
    % Constructor
    % -----------------------------------------------------------------
        
        %% ProcessedData
        function obj = ProcessedData( DeviceFileIn, IngestParamsIn )
        %#ProcessedData is the constructor for the object
        %#
        %# SYNOPSIS obj = ProcessedData( DeviceFileIn, IngestParamsIn )
        %# INPUT  DeviceFileIn   - a ACDeviceFile object
        %#        IngestParamsIn - the params struct from readIngestParameters
        %# OUTPUT obj            - the object
        %#
        
            %%% Pre-initialization %%%
            % Any code not using output argument (obj)
            
            %%% Object Initialization %%%
            % Call superclass constructor before accessing object
            % This statement cannot be conditionalized
            
            %%% Post-initialization %%%
            % Any code, including access to the object
            obj.L = log4m.getLogger();
            obj.L.debug('ProcessedData.ProcessedData()','Created object');
            
            if nargin > 0 
                % set device file
                if isa( DeviceFileIn, 'ACDeviceFile' )
                    obj.meta.DeviceFile = DeviceFileIn;
                else
                    error('Supply DeviceFile object')
                end
                
                if isa( IngestParamsIn, 'struct' )
                    obj.meta.Params = IngestParamsIn;
                else
                    error('Supply Ingest Params')
                end
            end   % if nargin > 0
            
            % set number of wavelengths from device file
            obj.meta.numWavelengths = str2num(obj.meta.DeviceFile.NumberWavelengths);
       
            keySet = {'raw', ...   %L1
                'preprocessed', ...%L2
                'binned', ...      %L3
                'filtered', ...    %L4
                'particulate', ... %L5
                'matchedWL', ...   %L6
                'corrected', ...   %L7
                'unsmoothed', ...  %L8
                'flagsApplied'}; %L9 MAYBE

            valueSet = {'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8', 'L9'};
            
            obj.levelsMap = containers.Map(keySet, valueSet);

        end  %#ProcessedData
        
    %% ----------------------------------------------------------------
    % Getter and Setter for Variables
    % -----------------------------------------------------------------
        %% setVar
        function obj = setVar( obj, varNameIn, levelNameIn, dataNameIn, dataIn )
        %#setVar sets any name/data pair into the right level in .var
        %#
        %# SYNOPSIS obj = setVar( obj, varNameIn, levelNameIn, dataNameIn, dataIn )
        %# INPUT  obj            - the object
        %#        varNameIn      - the var of the name/data pair being set, i.e. "ap"
        %#        levelNameIn    - the name of the level to set the data in, i.e. "corrected"
        %#        dataNameIn     - the name of the actual data, i.e. "timestamps"
        %#        dataIn         - the actual data
        %# OUTPUT obj            - the object
        %#
            
% debug statements in this code are printing regardless?            
%             obj.L.debug('ProcessedData.setVar()', 'in method');

            level = obj.levelsMap(levelNameIn);
            obj.var.(varNameIn).(level).(dataNameIn) = dataIn;
            
            %find index number for this level
            levelInt = level(:,2);  % get just the number off
            levelInt = str2num(levelInt);
           

            if isempty(obj.levelsFlags)
                
                % if no flags exist at all yet -- create flags and index
%                 obj.L.debug('ProcessedData.setVar()', 'have no flags at all');
                levelFlags = 1:9;
                levelIndex = zeros(size(levelFlags));
                levelIndex = logical(levelIndex);
                obj.levelsFlags.(varNameIn).levelFlags = levelFlags;
                obj.levelsFlags.(varNameIn).levelIndex = levelIndex;
                obj.levelsFlags.(varNameIn).levelIndex(levelInt) = true;
                
            else
%                 obj.L.debug('ProcessedData.setVar()','we do have flags');
                
                if ~isfield(obj.levelsFlags, varNameIn )
                    
                    % it's not empty - but need to check if we have this level
%                     obj.L.debug('ProcessedData.setVar()','have flags but not for this data');
                    levelFlags = 1:9;
                    levelIndex = zeros(size(levelFlags));
                    levelIndex = logical(levelIndex);
                    obj.levelsFlags.(varNameIn).levelFlags = levelFlags;
                    obj.levelsFlags.(varNameIn).levelIndex = levelIndex;
                    obj.levelsFlags.(varNameIn).levelIndex(levelInt) = true;
                    
                else
                    
                    % just update correct level
%                     obj.L.debug('ProcessedData.setVar()','have flags for this data');
                    obj.levelsFlags.(varNameIn).levelIndex(levelInt) = true;
                end;
            end
                
        end; %#setVar
        
        %% getVar
        function [varOut] = getVar( obj, varargin )
        %#getVar gets any data out of the right level in .var
        %#
        %# SYNOPSIS [varOut] = getVar( obj, varargin )
        %# INPUT  obj            - the object
        %#        varNameIn      - the var of the name/data pair being set, i.e. "ap"
        %#        levelNameIn    - the name of the level to set the data in, i.e. "corrected"
        %#        dataNameIn     - the name of the actual data, i.e. "timestamps"
        %# OUTPUT varOut         - the data
        %#  
        
% debug statements in this code are printing regardless?               
%             obj.L.debug('ProcessedData.getVar()', 'in method');
            
            % check inputs
            % check varName is one that exists. Create if it doesn't exist?
            varNameIn = '';
            levelNameIn = '';
            dataNameIn = '';
            
            % loop through name/value pairs
            if (~isempty(varargin))
                iArg = 1;
                while iArg < nargin
                    if strcmpi(varargin{iArg}, 'name')
                        varNameIn = varargin{iArg+1};
                    elseif strcmpi(varargin{iArg}, 'level')
                        levelNameIn = varargin{iArg+1};
                    elseif strcmpi(varargin{iArg}, 'data')
                        dataNameIn = varargin{iArg+1};
                    else
                        obj.L.error('ProcessedData.getVar', 'invalid argument');
                    end;
                    iArg = iArg + 2;
                end;  % while loop
            else
                obj.L.error('ProcessedData.getVar', 'no argument');
            end;
            
            % get variable - by name, level and datatype
            if ~isempty(varNameIn) && ~isempty(levelNameIn) && ~isempty(dataNameIn)
                
                % lookup level
                level = obj.levelsMap(levelNameIn);
                
                % check this data field exists for this level of data
                if isfield(obj.var.(varNameIn).(level), dataNameIn )
                    varOut = obj.var.(varNameIn).(level).(dataNameIn);
%                     obj.L.debug('ProcessedData.getVar()', ...
%                         sprintf('level: %s', level));
                else
                    obj.L.error('ProcessedData.getVar', 'invalid data name');
                end;
                
            elseif ~isempty(varNameIn) && ~isempty(levelNameIn) && isempty(dataNameIn)
%                 obj.L.debug('ProcessedData.getVar()', ...
%                     'have both name and level, don''t need data -- getting specific level');
                
                % lookup level
                level = obj.levelsMap(levelNameIn);

                if isfield(obj.var.(varNameIn), level )
                    varOut = obj.var.(varNameIn).(level);
                else
                    obj.L.error('ProcessedData.getVar', 'invalid level');
                end;

            elseif ~isempty(varNameIn) && isempty(levelNameIn)  && isempty(dataNameIn)
%                 obj.L.debug('ProcessedData.getVar()', ...
%                     'have name, don''t need data -- need to find most recent level');
                
                % find most recent level
                idx = obj.levelsFlags.(varNameIn).levelIndex;
                thisLevels = obj.levelsFlags.(varNameIn).levelFlags(idx);
                
                % if we have data for this var:
                if ~isempty(thisLevels)
                    maxLevel = max(thisLevels);
                    level = sprintf('L%u', maxLevel);
                    varOut = obj.var.(varNameIn).(level);
                else
                    obj.L.error('ProcessedData.getVar()', 'no data for this variable');
                end;
                

            elseif ~isempty(varNameIn) && isempty(levelNameIn)  && ~isempty(dataNameIn)
            % if we have a variable name and a data field name, but no
            % specific level, get the highest level we have for this
            % data field for this variable                
                
                idx = obj.levelsFlags.(varNameIn).levelIndex;
                thisLevels = obj.levelsFlags.(varNameIn).levelFlags(idx);
                
                % if we have data for this var:
                if ~isempty(thisLevels)
                    minLevel = min(thisLevels);
                    maxLevel = max(thisLevels);
                    for iLevel = minLevel:maxLevel
                        level = sprintf('L%u', iLevel);
                        levelExists = isfield(obj.var.(varNameIn), level);
                        if levelExists 
                            dataExists = isfield(obj.var.(varNameIn).(level), dataNameIn );
                            if dataExists % at this level
                                varOut = obj.var.(varNameIn).(level).(dataNameIn);
%                                 obj.L.debug('ProcessedData.getVar()', ...
%                                     sprintf('setting data to level: %s', level));

                            else
%                                 obj.L.debug('ProcessedData.getVar()', ...
%                                     'data doesn''t exist at this level');
                            end;
                        else
%                             obj.L.debug('ProcessedData.getVar()', ...
%                                 'level doesn''t exist');
                        end;
                    end;   % for loop through levels
                else
                    obj.L.error('ProcessedData.getVar()', 'no data for this variable');
                end;  %isempty(thisLevels)              
            else
                obj.L.error('ProcessedData.getVar', 'problem');
            end; %~isempty(varNameIn) && ~isempty(levelNameIn) && ~isempty(dataNameIn)
                
        end; %#getVar
      
    %% --------------------------------------------------------------------
    % Processing methods
    % ---------------------------------------------------------------------

        %% processBins
        function obj = processBins(obj, varargin)
        %#processBins - Take L2 "preprocessed" data and bin it -> L3 "binned".  
        %#If no argsin are provided, it will process both a/c and TSW/FSW
        %#            - calls processFSW()
        %# SYNOPSIS obj = processBins(obj, varargin)
        %# INPUT  obj            - the object
        %#        dataType1      - the primary type of the data being binned, i.e. "a" or "c"
        %#        dataType2      - the secondary type of the data being binned, i.e. TSW or FSW
        %# OUTPUT obj            - the object
        %# 

            obj.L.info('ProcessBins','Start of Method');
           
            % check varargin
            dataType1 = '';
            dataType2 = '';

            % parse optional input parameters
            if (~isempty(varargin))
                for iArgs = 1:length(varargin)
                    switch varargin{iArgs}
                        case {'a'}
                            dataType1 = {'a'};
                        case {'c'}
                            dataType1 = {'c'};
                        case {'TSW'}
                            dataType2 = {'TSW'};
                        case {'FSW'}
                            dataType2 = {'FSW'};
                        otherwise
                            
                           obj.L.error('ProcessedData.processBins', 'Invalid argument');
                    end   % switch
                end   % for
            else
                obj.L.debug('ProcessedData.processBins', 'no args in');
                dataType1 = {'a','c'};
                dataType2 = {'TSW','FSW'};
            end;   % if varargin is empty
            
            for iType1 = 1:length(dataType1)
                  for iType2 = 1:length(dataType2)

                      obj.L.info('ProcessedData.processBins', sprintf('type: %s', dataType1{iType1}));
                      obj.L.info('ProcessedData.processBins', sprintf('type: %s', dataType2{iType2}));
                      
                      thisType = sprintf('%s%s', dataType1{iType1}, dataType2{iType2});
            
                      [binnedTime, numberBins, binIndexNumbers, all_sample_size, ....
                          all_mean, all_median, all_std, all_variance, ...
                          all_variability, all_variability2,...
                          sample_size, ...
                          mean, median, std, variance, variability, variability2] = ...
                          processFSWTSW( obj, obj.getVar('name', thisType), ...
                          obj.meta.Params.PROCESS.BIN_SIZE );

                      if all(isnan(std))
                          obj.L.debug('ProcessBins', 'all std nan');
                      else
                          obj.L.debug('ProcessBins', 'all std not nan');
                      end;
                           
                    obj.setVar( thisType, 'binned', 'binnedTime', binnedTime);
                    obj.setVar( thisType, 'binned', 'numberBins', numberBins);
                    obj.setVar( thisType, 'binned', 'binIndexNumbers', binIndexNumbers);

                    
                    % these are the 'despiked' data
                    obj.setVar( thisType, 'binned', 'sample_size', sample_size);
                    obj.setVar( thisType, 'binned', 'mean', mean);
                    obj.setVar( thisType, 'binned', 'median', median);                    
                    obj.setVar( thisType, 'binned', 'std', std);
                    obj.setVar( thisType, 'binned', 'variance', variance);                    
                    obj.setVar( thisType, 'binned', 'variability', variability);
                    obj.setVar( thisType, 'binned', 'variability2', variability2);
                    
                    obj.setVar( thisType, 'binned', 'all_sample_size', all_sample_size);
                    obj.setVar( thisType, 'binned', 'all_mean', all_mean);
                    obj.setVar( thisType, 'binned', 'all_median', all_median);                    
                    obj.setVar( thisType, 'binned', 'all_std', all_std);
                    obj.setVar( thisType, 'binned', 'all_variance', all_variance);                    
                    obj.setVar( thisType, 'binned', 'all_variability', all_variability);
                    obj.setVar( thisType, 'binned', 'all_variability2', all_variability2);
                        
                  end   %for iType2
            end   %for iType1
            
            obj.L.info('ProcessedData.ProcessBins:','End of Method');         
        end  %#processBins

        %% processFSWTSW 
        function  [binned_time, numberBins, binIndexNumbers, ...
                all_sample_size, ...
                all_mean, all_median, all_std, all_variance, ...
                all_variability, all_variability2, ...
                despiked_sample_size, ...
                despiked_mean, despiked_median, despiked_std, ...
                despiked_variance, ...
                despiked_variability, despiked_variability2] ...
                = processFSWTSW(obj, dataIn, binSizeIn )
        %#processFSW: actually processes BOTH FSW and TSW for a AND c  
        %#If no argsin are provided, it will process both a/c and TSW/FSW
        %#
        %# SYNOPSIS obj = processBins(obj, varargin)
        %# INPUT  obj                  - the object
        %#        dataIn               - the actual data to bin
        %#        binSizeIn            - a datenum bin size to use, 
        %#                               i.e. datenum(0,0,0,0,1,0) for 1 minute
        %# OUTPUT 
        %#        binned_time          - the binned timestamps
        %#        numberBins           - the number of bins created
        %#        binIndexNumbers      - a vector of index numbers for the bins for the raw data
        
            
            
            
            obj.L.info('ProcessedData.processFSW','Start of Method');  

            data = dataIn.data;
            time = dataIn.timestamps;
            obj.L.debug('ProcessedData.processFSW', sprintf('First timestamp to be binned: %s', datestr(time(1))));
            obj.L.debug('ProcessedData.processFSW', sprintf('Last  timestamp to be binned: %s', datestr(time(end))));
            
            binSize = binSizeIn;
            UPPER_PERCENTILE = obj.meta.Params.PROCESS.UPPER_PERCENTILE;
            LOWER_PERCENTILE = obj.meta.Params.PROCESS.LOWER_PERCENTILE;
            
            % CHANGED FOR NAAMES -- MISMATCH BETWEEN DEVICE FILE AND DATA
            [rows,cols] = size(data);
            obj.L.debug('ProcessedData.processFSW',sprintf('size of data: %u x %u', rows, cols));
            numberWavelengths = cols;
            obj.L.debug('ProcessedData.processFSW',sprintf('numberWavelengths: %u', numberWavelengths));
            % END CHANGE FOR NAAMES
            
            % Binned time stamps
            % convert time datenums to datevec type
            time_vec = datevec(time);
            % create a datevec of the first time stamp
            start_time_vec = time_vec(1,:);
            % set the seconds to zero (just have minute)
            start_time_vec(6) = 0;
            % convert it back to a datenum
            start_time = datenum(start_time_vec);
            % create a datevec of the last timestamp
            end_time_vec = time_vec(end,:);
            % set the seconds to zero (just have minutes)
            end_time_vec(6) = 0;
            % convert it back to da datenum
            end_time = datenum(end_time_vec);
            
            % create bins of minutes, starting at start time, ending at end time, with
            % one minute intervals
            binned_time = start_time:binSize:end_time;
            numberBins = numel(binned_time);
            
            obj.L.debug('ProcessedData.processFSW', sprintf('Number of bins: %u', numberBins));
            obj.L.debug('ProcessedData.processFSW', sprintf('First bin: %s', datestr(binned_time(1))));
            obj.L.debug('ProcessedData.processFSW', sprintf('Last  bin: %s', datestr(binned_time(end))));

            % create an index of the bins
            binIndexNumbers = zeros(size(time));
            binIndexNumbers(:,:) = NaN;

            % sample_size = number of timestamps in Bin
            all_sample_size = zeros(numberBins, numberWavelengths);
            all_sample_size(:,:) = NaN;
            
            % set up matrices for statistics:
            % one set for all data, one set for despiked
            all_mean = zeros(numberBins, numberWavelengths);
            all_mean(:,:) = NaN;
            
            all_median = zeros(numberBins, numberWavelengths);
            all_median(:,:) = NaN;
             
            all_std = zeros(numberBins, numberWavelengths);
            all_std(:,:) = NaN;
            
            all_variance = zeros(numberBins, numberWavelengths);
            all_variance(:,:) = NaN;
            
            all_variability = zeros(numberBins, numberWavelengths);
            all_variability(:,:) = NaN;
            
            all_variability2 = zeros(numberBins, numberWavelengths);
            all_variability2(:,:) = NaN;
            
            % sample_size = number of timestamps in Bin
            despiked_sample_size = zeros(numberBins, numberWavelengths);
            despiked_sample_size(:,:) = NaN;     
            
            despiked_mean = zeros(numberBins, numberWavelengths);
            despiked_mean(:,:) = NaN;
            
            despiked_median = zeros(numberBins, numberWavelengths);
            despiked_median(:,:) = NaN;
             
            despiked_std = zeros(numberBins, numberWavelengths);
            despiked_std(:,:) = NaN;
            
            despiked_variance = zeros(numberBins, numberWavelengths);
            despiked_variance(:,:) = NaN;
            
            despiked_variability = zeros(numberBins, numberWavelengths);
            despiked_variability(:,:) = NaN;
            
            despiked_variability2 = zeros(numberBins, numberWavelengths);
            despiked_variability2(:,:) = NaN;

            % loop through each bin (currently minute)
            for iBin=1:numel(binned_time)

                % get an index of which timestamps belong in this bin
                if iBin < numel(binned_time)
                    % its any bin other than the last one
                    thisBinTimestampIndex = ( time >= binned_time(iBin) ) & ( time < binned_time(iBin+1) );
                else
                    % its the last timestamp
                    thisBinTimestampIndex = ( time >= binned_time(iBin) );
                end

                % Fill in the index number of this bin for future processing
                binIndexNumbers(thisBinTimestampIndex) = iBin;

                % Gather the data for this bin
                thisBinData = data(thisBinTimestampIndex,:);
                [numRows, numCols] = size(thisBinData);

%                 all_sample_size(iBin,:) = numRows;
                %sample size is all the rows (,2) in this matrix that don't
                %have all nans
                all_sample_size(iBin,:) = sum(~all(isnan(thisBinData),2));
                
                % moved up here 6/7/16
                all_mean(iBin,:) = nanmean(thisBinData, 1);
                all_median(iBin,:) = nanmedian(thisBinData, 1);
                all_std(iBin,:) = nanstd(thisBinData, 1);
                all_variance(iBin,:) = nanvar(thisBinData, 1);
                
                % LOWER_PERCENTILE is 2.5% default 
                % UPPER_PERCENTILE is 97.5% default
                all_variability(iBin,:) = ...
                     ( prctile(thisBinData, UPPER_PERCENTILE) - prctile(thisBinData,LOWER_PERCENTILE) )/2;
                all_variability2(iBin,:) = ...
                     ( prctile(thisBinData,84) - prctile(thisBinData,16) )/2;

                % check not ALL of these values are NaN
                if ~all(all(isnan(thisBinData)))
                    
                    % 1.  FIND PERCENTILES
                    %     Find data in between 2.5% and 97.5% in this bin
                    %     prctile will return a matrix - one column per column of data 
                    %     first row for LOWER_PERCENTILE value
                    %     second row for UPPER_PERCENTILE value
                    thisBinPercentiles = prctile( thisBinData, [LOWER_PERCENTILE UPPER_PERCENTILE], 1);

                    % set up an index of good data to use
                    thisBinPercCheckIndex =  zeros( numRows, numCols);
                    thisBinPercCheckIndex(:,:) = NaN;

                    % for each column of data (for each wavelength)
                    %     create an index of data to use in further calculations
                    %     if it passes the check, assign a 1 to the index
                    for iCol = 1:numCols

                        lowerPercentile = thisBinPercentiles(1, iCol);
                        upperPercentile = thisBinPercentiles(2, iCol);

                        thisBinPercCheckIndex(:,iCol) = ( thisBinData(:,iCol) >=  lowerPercentile) & ...
                            ( thisBinData(:, iCol) <= upperPercentile );

                    end;

                    thisBinPercCheckIndex = logical(thisBinPercCheckIndex);
                    
                    %create a matrix, same dimensions as data, to copy data with index
                    %applies into.
                    thisBinPercCheckGood = thisBinData;
                    thisBinPercCheckGood(:,:) = NaN;

                    % for each location, copy in the good data
                    thisBinPercCheckGood(thisBinPercCheckIndex) = thisBinData(thisBinPercCheckIndex);
                     
                    despiked_sample_size(iBin,:) = sum(~all(isnan(thisBinPercCheckGood),2));
                    
                    % 2.  STATISTICS
                    %     use nan(STAT), i.e. nanmean:  "For matrices X, nanmean(X) is a row vector of
                    %     column means, once NaN values are removed."
                    %     thisBinData is N number of timestamps, with 83 columns.
                    %     We want column means -- one for each WL
                    %     Assign to row "i" of bin_data_mean
                    
                    % Only use data that met percentile criteria above
                    % (thisBinPercCheckGood)

                    despiked_mean(iBin,:) = nanmean(thisBinPercCheckGood, 1);
                    despiked_median(iBin,:) = nanmedian(thisBinPercCheckGood, 1);                    
                    despiked_std(iBin,:) = nanstd(thisBinPercCheckGood, 1);
                    despiked_variance(iBin,:) = nanvar(thisBinPercCheckGood,1);
 
                    thisBinFinalPercentiles = prctile( thisBinPercCheckGood, [LOWER_PERCENTILE UPPER_PERCENTILE], 1);

                    for iCol = 1:numCols
                        lowerPercentile = thisBinFinalPercentiles(1, iCol);
                        upperPercentile = thisBinFinalPercentiles(2, iCol);
                        % (97.%-2.5%)/2           
                        despiked_variability(iBin, iCol) = (upperPercentile-lowerPercentile)/2;
                        
                    end;
                    despiked_variability2(iBin,:) = ...
                        ( prctile(thisBinPercCheckGood,84) - prctile(thisBinPercCheckGood,16) )/2;
                    
                   
                else  %  ~all(all(isnan(thisBinData)))
                     obj.L.debug('ACData.processFSW', sprintf('Bin %u has no valid data', iBin));
                end   % end check for valid data

            end;    %for each bin
            obj.L.info('processFSW','End of Method');           
        end   %#processFSW
        
        %% findFSWBinMedians
        function obj = findFSWBinMedians(obj, varargin)
        %#findFSWBinMedians - Interpolate Filtered Data for a and c: L3->L4
        %#                  - calls obj.findBinMeds
        %#
        %# SYNOPSIS obj = findFSWBinMedians(obj, varargin)
        %# INPUT obj       : the object
        %#       dataType1 : the primary type of the data being binned, i.e. "a" or "c"
        %# OUTPUT obj: the object
          
            % check varargin
            dataType1 = '';

            % parse optional input parameters
            if (~isempty(varargin))
                for iArgs = 1:length(varargin)
                    switch varargin{iArgs}
                        case {'a'}
                            dataType1 = {'a'};
                        case {'c'}
                            dataType1 = {'c'};
                        otherwise
                           obj.L.error('ProcessedData.findFSWBinMedian', 'Invalid argument');
                    end   % switch
                end   % for
            else
                obj.L.debug('ProcessedData.findFSWBinMedian', 'no args in');
                dataType1 = {'a','c'};
            end;   % if varargin is empty   
            
            for iType1 = 1:length(dataType1)  
                obj.L.debug('ProcessedData.findFSWBinMedian', sprintf('type: %s', dataType1{iType1}));                
                thisType = sprintf('%sFSW', dataType1{iType1});   
                binMethod = obj.meta.Params.PROCESS.STAT_FOR_BINNING;
                
                % get most recent data of the type given
                filtered_bins = obj.getVar('name', thisType, 'level', 'binned', 'data', binMethod);   %1502x83

                obj.L.debug('ProcessedData.findFSWBinMedian',...
                  sprintf('size filtered bins: %u x %u', size(filtered_bins)));

                % run method
                [binMedians] = obj.findBinMeds(filtered_bins);                
                %set variables
                obj.setVar( thisType, 'filtered', 'interpolatedSectionMedians', binMedians);

            end; %# for iType1 = 1:length(dataType1)
            obj.L.info('ProcessedData.findFSWBinMedian','End of Method');
        end; % #findFSWBinMedian(obj, varargin)            
        
        %% findBinMeds        
        function [binMediansOut] = findBinMeds(obj, filteredBinsIn)
        %#findBinMeds - find the median of each section of FSW bins, 
        %#and assign the median to the middle timestamp
        %#
        %# SYNOPSIS [binMediansOut] = findBinMeds(obj, filteredBinsIn)
        %# INPUT obj: the object
        %#       filteredBinsIn: FSW data that has been binned
        %# OUTPUT binMediansOut: data including one bin for each section of 
        %#                       bins, where value = median of that section
        %# 
            % create index to bins which are all NaNs: i.e. the TSW bins
            nanDataBinIndex = all(isnan(filteredBinsIn),2);  %1502x83

            % create an index of which bins have filtered data
            filtBinIndex = ~nanDataBinIndex;
            
            % calc size of data
            [numBins,numCols] = size(filteredBinsIn); 
             
            % create new matrix of medians, same size as original data
            interpolatedSectionMedians = zeros(numBins, numCols);
            interpolatedSectionMedians(:,:) = NaN;   
             
            % find where sections of FSW/TSW switch:
            offsetfiltBinIndex = zeros(size(filtBinIndex));
            offsetfiltBinIndex(1,:) = NaN;
            offsetfiltBinIndex(2:end,:) = filtBinIndex(1:end-1);
            transitionIndex = filtBinIndex - offsetfiltBinIndex;
            startSectionIndex  = find(transitionIndex == 1);
            endSectionIndex = find(transitionIndex == -1);
            % decrease section by one
            endSectionIndex = endSectionIndex -1;
            numSections = size(startSectionIndex);
            
            if length(startSectionIndex) > length(endSectionIndex)
                %ending in middle of FSW
                newEndSectionIndex = zeros(length(endSectionIndex)+1,1);
                newEndSectionIndex(1:end-1) = endSectionIndex(:,:);
                newEndSectionIndex(end) = length(filtBinIndex);
                endSectionIndex = newEndSectionIndex;
            end;

            % for each section
            for iSection = 1:numSections;

                 % create a temporary section the size of the current
                 % section of interpolated data
                 endSectionIndex(iSection);
                 startSectionIndex(iSection);
                 currSection = zeros( endSectionIndex(iSection) - startSectionIndex(iSection) + 1, numCols);
                 currSection(:,:) = NaN;
                 
                 % copy in the data from the original data
                 istart = startSectionIndex(iSection);
                 iend = endSectionIndex(iSection);
                 
                % find size of currSection                 
                 currSection(:,:) = filteredBinsIn( istart:iend, :) ;
                
                [numBinsInSection,~] = size(currSection);
                halfwayBin = startSectionIndex(iSection) + floor(numBinsInSection/2);
                
                obj.L.debug('ProcessedData.calculateInterpolatedBinUncertainty',...
                 sprintf('Size: %u; Halfway bin: %u', numBinsInSection, halfwayBin));

                obj.L.debug('ProcessedData.calculateInterpolatedBinUncertainty',...
                 sprintf('Section Number: %u; start: %u; end %u', ...
                 iSection, startSectionIndex(iSection), endSectionIndex(iSection)));

                medianCurrSection = nanmedian(currSection(:,:),1);
                 obj.L.debug('ProcessedData.calculateInterpolatedBinUncertainty',...
                     sprintf('Section Meidan: %u; median size: %u x %u',medianCurrSection, size(medianCurrSection)));                                              
                
                
                 for iBin = startSectionIndex(iSection):endSectionIndex(iSection)
                     for iCol = 1:numCols                 
                         currMax = max(currSection(:,iCol));
                         currMin = min(currSection(:,iCol));
                         % check interpolated bins
                        if abs(currMax-currMin) > 0.005
                         obj.L.error('ProcessedData.calculateInterpolatedBinUncertainty',...
                             sprintf('max-min > 0.005 = %u -- HIGH VAR in FSW Sec %u -- MAKE INTERVAL SMALLER?', ...
                             (abs(currMax-currMin) > 0.005), iBin));
                        end;
                                
                         % copy in median to middle bin
                         if iBin == halfwayBin
                             interpolatedSectionMedians(iBin,iCol) = medianCurrSection(:,iCol);
                         end  
                         
                     end %# for iCol = 1:numCols 
                 end %# for iBin = startSectionIndex(iSection):endSectionIndex(iSection)
            end
             
             binMediansOut = interpolatedSectionMedians;
             
         end %# findFilteredBinMedian

        %% interpolateFiltered
        function obj = interpolateFiltered(obj, varargin)
        %#interpolateFiltered - Interpolate Filtered Data for a and c: L3->L4
        %#                    - calls obj.interpolateFSW
        %#                    - calls obj.calculateInterpolatedBinUncertainty
        %#
        %# SYNOPSIS obj = interpolateFiltered(obj, varargin)
        %# INPUT obj: the object
        %#        dataType1      - the primary type of the data being binned, i.e. "a" or "c"
        %#        dataType2      - the secondary type of the data being binned, i.e. TSW or FSW        
        %# OUTPUT obj: the object
        %#            
        
            obj.L.info('ProcessedData.interpolateFiltered','Start of Method');
            
            % check varargin
            dataType1 = '';
            % only interpolate filtered data
            dataType2 = {'FSW'};

            % parse optional input parameters
            if (~isempty(varargin))
                for iArgs = 1:length(varargin)
                    switch varargin{iArgs}
                        case {'a'}
                            dataType1 = {'a'};
                        case {'c'}
                            dataType1 = {'c'};
                        otherwise
                            
                           obj.L.error('ProcessedData.interpolateFiltered', 'Invalid argument');
                    end   % switch
                end   % for
            else
                obj.L.debug('ProcessedData.interpolateFiltered', 'no args in');
                dataType1 = {'a','c'};
            end;   % if varargin is empty
            
            for iType1 = 1:length(dataType1)
                  for iType2 = 1:length(dataType2)

                      obj.L.debug('ProcessedData.interpolateFiltered', sprintf('type: %s', dataType1{iType1}));
                      obj.L.debug('ProcessedData.interpolateFiltered', sprintf('type: %s', dataType2{iType2}));
                      
                      thisType = sprintf('%s%s', dataType1{iType1}, dataType2{iType2});
                      
                     % get most recent data of the type given
                      filtered_bins = obj.getVar('name', thisType, 'level', ...
                          'filtered', 'data', 'interpolatedSectionMedians');   %1502x83
                      time_bins = obj.getVar('name', thisType, 'data', 'binnedTime'); %1x1502  
                      
                      obj.L.debug('ProcessedData.interpolateFiltered',...
                          sprintf('size filtered bins: %u x %u', size(filtered_bins)));
                      obj.L.debug('ProcessedData.interpolateFiltered',...                      
                          sprintf('size time bins: %u x %u', size(time_bins)));
                      
                      time_bins = time_bins';
                      
                    % run method
                    [interpolatedFSWData, interpolatedBinIndex] = ...
                        obj.interpolateFSW(time_bins, filtered_bins, 'extrap');
                      
                    %set variables
                    obj.setVar( thisType, 'filtered', 'interpolatedData', interpolatedFSWData);
                    obj.setVar( thisType, 'filtered', 'interpolatedBinIndex', interpolatedBinIndex);
                    obj.setVar( thisType, 'filtered', 'binnedTime', time_bins);
                    
                    % calculate bin variability
                    [uncertainty, ~ ] = obj.calculateInterpolatedBinUncertainty( ...
                        interpolatedFSWData, interpolatedBinIndex);
                    
                    % set variability
                    obj.setVar( thisType, 'filtered', 'interpolatedUncertainty', uncertainty);
                
                  end;
            end;
            
            obj.L.info('ProcessedData.interpolateFiltered','End of Method');
        end   %#interpolateFiltered
                      
        %% calcTSWUncertainty
        function obj = calcTSWUncertainty(obj, varargin)
            obj.L.info('ProcessedData.calcTSWUncertainty','Start of Method');

            % check varargin
            dataType1 = '';
            % only interpolate filtered data
            dataType2 = {'TSW'};

            % parse optional input parameters
            if (~isempty(varargin))
                for iArgs = 1:length(varargin)
                    switch varargin{iArgs}
                        case {'a'}
                            dataType1 = {'a'};
                        case {'c'}
                            dataType1 = {'c'};
                        otherwise
                         
                           obj.L.error('ProcessedData.calcTSWUncertainty', 'Invalid argument');
                    end   % switch
                end   % for
            else
                obj.L.debug('ProcessedData.calcTSWUncertainty', 'no args in');
                dataType1 = {'a','c'};
            end;   % if varargin is empty            
            
            for iType1 = 1:length(dataType1)
                  for iType2 = 1:length(dataType2)

                      obj.L.debug('ProcessedData.calcTSWUncertainty', sprintf('type: %s', dataType1{iType1}));
                      obj.L.debug('ProcessedData.calcTSWUncertainty', sprintf('type: %s', dataType2{iType2}));
                      
                      thisType = sprintf('%s%s', dataType1{iType1}, dataType2{iType2});
                      
                       % get most recent data of the type given
                      TSW_STD = obj.getVar('name', thisType, 'level', ...
                          'binned', 'data', 'std');   %1502x83
                      TSW_Var = obj.getVar('name', thisType, 'level', ...
                          'binned', 'data', 'variability');   %1502x83
                      
                      if strcmpi(obj.meta.Params.PROCESS.BIN_METHOD, 'mean')
                           binUncertainty = TSW_STD./sqrt(obj.meta.Params.PROCESS.UNCERTAINTY_N);
                      elseif strcmpi(obj.meta.Params.PROCESS.BIN_METHOD, 'median')
                           binUncertainty = TSW_Var./sqrt(obj.meta.Params.PROCESS.UNCERTAINTY_N);
                      else
                           obj.L.error('ProcessingManager', 'unexpected bin method');
                      end;
                      
                      obj.setVar( thisType, 'binned', 'uncertainty', binUncertainty);
                      
                  end;  %for iType2
            end; %for iType1
        end;  %calcTSWUncertainty
       
        %% interpolateFSW
        function [interpolatedFSWData, interpolatedBinIndex] = interpolateFSW(obj, binTimestampsIn, binnedFSWDataIn, extrapIn)
        %#interpolateFSW - use INTERP1 to create interpolated data between FSW
        %#
        %# SYNOPSIS [interpolatedFSWData, interpolatedBinIndex] = interpolateFSW(obj, binTimestampsIn, binnedFSWDataIn)
        %# INPUT obj: the object
        %#       binTimestampsIn
        %#       binnedFSWDataIn
        %# OUTPUT interpolatedFSWData
        %#        interpolatedBinIndex
        %#            
        
            obj.L.info('ProcessedData.interpolateFSW','Start of Method');
            
            % STEP 1:
            % When interpolating, the array that interp1 uses to find the values can't
            % have any NaNs so create a copy of your filtered array to fill in the blanks.

            % make copy of binned data passed in, NaNs will be blanked out of this
            FSWBinsNoNans = binnedFSWDataIn;  %1502x83
            
            % make copy of binned time passed in, NaNs will be blanked out of this 
            binTimestampsNoNans = binTimestampsIn;
            
            % make another copy of binned timestamps, to use for
            % interpolation
            binTimestamps = binTimestampsIn;
            
            % create index to bins which are all NaNs: i.e. the TSW bins
            nanDataBinIndex = all(isnan(FSWBinsNoNans),2);  %1502x83

            % what is this for?
            filtBinIndex = ~nanDataBinIndex;

            % blank any bins that are totally nan
            FSWBinsNoNans(nanDataBinIndex, :) = [];  % new array with just data
            binTimestampsNoNans(nanDataBinIndex, :) = [];  % corresponding time_stamps

            % STEP 2:
            % make copy of filtered_data for interpolated data to go INTO
            % this is copy of original so WILL contain NaNs.
            % new_filtered_data = filtered_bins;
            if strcmp(extrapIn, 'extrap')
                interpolatedFSWData = interp1(binTimestampsNoNans, FSWBinsNoNans, binTimestamps, ...
                    'linear', 'extrap');   
            else
                interpolatedFSWData = interp1(binTimestampsNoNans, FSWBinsNoNans, binTimestamps, ...
                    'linear'); 
            end;
            newDataIndex = ~all(isnan(interpolatedFSWData),2);

            % an index to the interpolated bins will be the OPPOSITE of:
            % bins that were Nan before Interp AND bins that are nan after
            % Interpolation
            interpolatedBinIndex = newDataIndex & ~filtBinIndex;
            
            obj.L.info('ProcessedData.interpolateFSW','End of Method');

        end   %end interpolateFSW 
        
        %% calculateInterpolatedBinUncertainty
        function [interpolatedBinUncertaintyOut, interpolatedBinUncertainty2Out] ...
                = calculateInterpolatedBinUncertainty( ...
                obj, interpolatedFSWDataIn, interpolatedBinIndexIn ) 
             
        %#calculateInterpolatedBinUncertainty - calculates interpolated bin variability as max-min / 2
        %#
        %# SYNOPSIS [interpolatedBinUncertaintyOut, interpolatedBinUncertainty2Out ] =  ...
        %#        calculateInterpolatedBinUncertainty(obj, interpolatedFSWDataIn, interpolatedBinIndexIn ) 
        %# INPUT  obj: the object
        %#        interpolatedFSWDataIn -
        %#        interpolatedBinIndexIn - 
        %# OUTPUT interpolatedBinUncertaintyOut
        %#        interpolatedBinUncertainty2Out
        %#            
            obj.L.info('ProcessedData.calculateInterpolatedBinUncertainty','Start of Method');
            interpolatedFSWData = interpolatedFSWDataIn;
            interpBinIndex = interpolatedBinIndexIn;

             [~,numCols] = size(interpolatedFSWData); 
             [numBins] = length(interpolatedFSWData);
             
             interpolatedBinUncertainty = zeros(numBins, numCols);
             interpolatedBinUncertainty(:,:) = NaN;

             interpolatedBinUncertainty2 = zeros(numBins, numCols);
             interpolatedBinUncertainty2(:,:) = NaN;
          
             % calculate variability for each interpolated bin.
             % variability is max of bins in that section - min of bins in
             % that section/2
             offsetInterpBinIndex = zeros(size(interpBinIndex));
             offsetInterpBinIndex(1,:) = NaN;
             offsetInterpBinIndex(2:end,:) = interpBinIndex(1:end-1);
             transitionIndex = interpBinIndex - offsetInterpBinIndex;
             
             %switched
             % find where the transitionIndex is 1, this will be where to
             % start a new seciton.
             startSectionIndex  = find(transitionIndex == 1);
             
              %find where transitinIndex is -1, this will mark end of
             %section
             endSectionIndex = find(transitionIndex == -1);
             
             startStartsAtBeginning = (length(find(startSectionIndex == 1)) > 0);
             endStartsAtBeginning = (length(find(endSectionIndex == 1)) > 0);
             startEndsAtEnd = (length(find(startSectionIndex == numBins)) > 0);
             endEndsAtEnd = (length(find(endSectionIndex == numBins)) > 0);
             
             if ~(startStartsAtBeginning || endStartsAtBeginning || startEndsAtEnd || endEndsAtEnd)
                 % none of them start/end at start or end
                 % create a new index, one longer than existing
                 newStartSectionIndex = zeros((length(startSectionIndex) + 1), 1);
                 newStartSectionIndex(1,1) = 1;  % set first one to start
                 % copy rest in
                 newStartSectionIndex(2:end,1) = startSectionIndex;

                 newEndSectionIndex = zeros((length(startSectionIndex) + 1),1);
                 % create a new index to account for end, copy in existing
                 % index and subtract one to move backwards
                 newEndSectionIndex(1:end-1,1) = endSectionIndex-1;
                 newEndSectionIndex(end,1) = length(transitionIndex);
             else
                 % one of them starts/ends at start/end
                 if startStartsAtBeginning
                     newStartSectionIndex = startSectionIndex;
                 else
                     % create a new index, one longer than existing
                     newStartSectionIndex = zeros((length(startSectionIndex) + 1), 1);
                     newStartSectionIndex(1,1) = 1;  % set first one to start
                     newStartSectionIndex(2:end,1) = startSectionIndex;
                 end;
                 if endEndsAtEnd

                     % create a new index to account for end, copy in existing
                     % index and subtract one to move backwards
                     newEndSectionIndex = endSectionIndex-1;
                 else
                     newEndSectionIndex = zeros((length(endSectionIndex) + 1),1);
                     % create a new index to account for end, copy in existing
                     % index and subtract one to move backwards
                     newEndSectionIndex(1:end-1,1) = endSectionIndex-1;
                     newEndSectionIndex(end,1) = length(transitionIndex); 
                 end;

             end
             
             % copy over to previous names so code runs
             startSectionIndex = newStartSectionIndex;
             endSectionIndex = newEndSectionIndex;
             numSections = length(startSectionIndex);
             % for each section
             for iSection = 1:numSections;
                 
                 % create a temporary section the size of the current
                 % section of interpolated data
                 currSection = zeros( endSectionIndex(iSection) - startSectionIndex(iSection) + 1, numCols);
                 currSection(:,:) = NaN;
                 
                 % copy in the data from the original interp data
                 istart = startSectionIndex(iSection);
                 iend = endSectionIndex(iSection);
                 currSection(:,:) = interpolatedFSWData( istart:iend, :) ;
                 
                 for iBin = startSectionIndex(iSection):endSectionIndex(iSection)
                     for iCol = 1:numCols
                         currMax = max(currSection(:,iCol));
                         currMin = min(currSection(:,iCol));
                         currVar = (currMax - currMin)/2;
                         interpolatedBinUncertainty( iBin, iCol ) = currVar;
                         
                         % do same calculation, but for 1st and last
                         % measurements
                         currSecStart = currSection(1,iCol);
                         currSecEnd = currSection(end,iCol);
                         currSecVar = abs(currSecStart-currSecEnd)/2;
                         interpolatedBinUncertainty2(iBin,iCol) = currSecVar;
                         
                     end
                 end
             end
             
             interpolatedBinUncertaintyOut = interpolatedBinUncertainty;
             interpolatedBinUncertainty2Out = interpolatedBinUncertainty2;
             
             obj.L.info('ProcessedData.calculateInterpolatedBinUncertainty','End of Method');            
         end   %# calculateInterpolatedBinUncertainty        
        
        %% calcParticulate
        function obj = calcParticulate( obj, varargin )
        %#calcParticulate calculate cp and remove timestamps beyond current YEARDAY
        %#
        %# SYNOPSIS calcParticulate( obj, varargin )
        %# INPUT  obj       - the object
        %#        dataType1 - the primary type of the data being binned, i.e. "a" or "c"
        %#                  - if no type is given, both will be run
        %# OUTPUT obj       - the object
        %#          
            obj.L.info('ProcessedData.calcParticulate','Start of Method');    

            % could have only a or c

            dataType1 = '';

            % parse optional input parameters
            if (~isempty(varargin))
                for iArgs = 1:length(varargin)
                    switch varargin{iArgs}
                        case {'a'}
                            dataType1 = {'a'};
                        case {'c'}
                            dataType1 = {'c'};
                        otherwise
                           obj.L.error('ProcessedData.calcParticulate', 'Invalid argument');
                    end   % switch
                end   % for
            else
                obj.L.debug('ProcessedData.calcParticulate', 'no args in');
                
                % by default, run both data types a & c
                dataType1 = {'a','c'};
                
            end;   % if varargin is empty
            
            % for a,c
            for iType1 = 1:length(dataType1)
                
                % set the types of data that come before computing p
                thisType = sprintf('%sp', dataType1{iType1} );
                tswType = sprintf('%sTSW', dataType1{iType1} );
                fswType = sprintf('%sFSW', dataType1{iType1} );
                binMethod = obj.meta.Params.PROCESS.STAT_FOR_BINNING;
 
                % -------------------------------------------------------
                % get variables
                % -------------------------------------------------------
                
                % get TSW bins and Interpolated FSW bins for calculating p
                tswBin = obj.getVar('name', tswType, 'level', 'binned', 'data', binMethod);
                interpBin = obj.getVar('name', fswType, 'level', 'filtered', 'data', 'interpolatedData');
                
                % get variables for calculating uncertainty
                tswBinUncertainty = obj.getVar('name', tswType, 'level', 'binned', 'data', 'uncertainty' );
                interpBinUncertainty = obj.getVar('name', fswType, 'level', 'filtered', 'data', 'interpolatedUncertainty' );
                      
                % get timestamps to asssign to cp/ap
                timestamps = obj.getVar('name', fswType, 'level', 'filtered', 'data', 'binnedTime');
          
                
                % check matlab hasn't flipped arbitrarily
                [r,c] = size(timestamps);
                if r<c
                    timestamps = timestamps';
                end;
                % added 3/3/16
                this_std = obj.getVar('name', tswType, 'level', 'binned', 'data' , 'std');
                
                % added 4/5/16
                % get bin count
                binCount = obj.getVar('name', tswType, 'level', 'binned', 'data', 'sample_size');
                
                                
                TSW_all_mean = obj.getVar('name', tswType, 'level', 'binned', 'data', 'all_mean'); 
                TSW_all_var = obj.getVar('name', tswType, 'level', 'binned', 'data', 'all_variance'); 
                FSW_all_mean = obj.getVar('name', fswType, 'level', 'binned', 'data', 'all_mean');
                
                TSW_despiked_mean = obj.getVar('name', tswType, 'level', 'binned', 'data', 'mean'); 
                TSW_despiked_var = obj.getVar('name', tswType, 'level', 'binned', 'data', 'variance');
                
                FSW_despiked_mean = obj.getVar('name', fswType, 'level', 'binned', 'data', 'mean'); 
                
                % -------------------------------------------------------
                % make calculations
                % -------------------------------------------------------
                
                % 1. calculate particulate = total - filtered
                pUncorr = tswBin - interpBin;
                      
                % 2.  calculate uncertainty
                % UNCERTAINTY =  InterpolatedBinUncertainty PLUS
                % TSWBinUncertainty
                totalUncertainty = interpBinUncertainty + tswBinUncertainty;
                
                % 3. calculate var/mean for all data and despiked data
%                 TSW.init_var/[ TSW.init_mean – FSW.init_mean ]
%                 TSW.var/[ TSW.mean – FSW.mean ]

                % first interpolate the means
                [FSW_interp_mean, ~] = obj.interpolateFSW(timestamps, FSW_all_mean, 'no');
                [FSW_interp_despiked_mean, ~] = obj.interpolateFSW(timestamps, FSW_despiked_mean, 'no');

                all_var_over_mean = TSW_all_var./(TSW_all_mean - FSW_interp_mean);
                despiked_var_over_mean = TSW_despiked_var./(TSW_despiked_mean - FSW_interp_despiked_mean);
                
                % 4.  Now that cp is calculated, remove timestamps from
                % data set that are beyond this yearday
                
                % create an index of timestamps to remove
                removeIndex = timestamps > datenum(obj.meta.Params.INGEST.YEAR, ...
                    obj.meta.Params.INGEST.MONTH, obj.meta.Params.INGEST.DAY, 23, 59, 59) | ...
                    timestamps < datenum(obj.meta.Params.INGEST.YEAR, ...
                    obj.meta.Params.INGEST.MONTH, obj.meta.Params.INGEST.DAY);
                
                
                obj.L.debug('ProcessedData.calcParticulate', ...
                    sprintf('size of timestamps before: %u %u', size(timestamps)));
                obj.L.debug('ProcessedData.calcParticulate', ...
                    sprintf('LAST TIMESTAMP before: %s', datestr(timestamps(end))));
                obj.L.debug('ProcessedData.calcParticulate', ...
                    sprintf('size of pUncorr before: %u %u', size(pUncorr)));
                obj.L.debug('ProcessedData.calcParticulate', ...
                    sprintf('size of timestamps before: %u %u', size(totalUncertainty)));                
                obj.L.debug('ProcessedData.calcParticulate', ...
                    sprintf('size of std before: %u %u', size(this_std)));
                
                % remove timestamps from timestamps
                timestamps = removerows(timestamps, removeIndex);
                pUncorr = removerows(pUncorr, removeIndex);
                totalUncertainty = removerows(totalUncertainty, removeIndex);
                this_std = removerows(this_std, removeIndex);
                binCount = removerows(binCount, removeIndex);
                all_var_over_mean = removerows(all_var_over_mean, removeIndex);
                despiked_var_over_mean = removerows(despiked_var_over_mean, removeIndex);
                
                obj.L.debug('ProcessedData.calcParticulate', ...
                    sprintf('size of timestamps after: %u %u', size(timestamps)));
                obj.L.debug('ProcessedData.calcParticulate', ...
                    sprintf('LAST TIMESTAMP after: %s', datestr(timestamps(end))));                
                obj.L.debug('ProcessedData.calcParticulate', ...
                    sprintf('size of pUncorr after: %u %u', size(pUncorr)));
                obj.L.debug('ProcessedData.calcParticulate', ...
                    sprintf('size of timestamps after: %u %u', size(totalUncertainty)));   
                obj.L.debug('ProcessedData.calcParticulate', ...
                    sprintf('size of std after: %u %u', size(this_std)));  
                
                % -------------------------------------------------------
                % set variables
                % -------------------------------------------------------
                
                % set remove index
                obj.setVar( thisType, 'particulate', 'removeIndex', removeIndex);
                
                % set data
                obj.setVar( thisType, 'particulate', 'data', pUncorr);
                
                % set uncertainty
                obj.setVar( thisType, 'particulate', 'uncertainty', totalUncertainty);
                
                % set timestamps
                obj.setVar( thisType, 'particulate', 'timestamps', timestamps);
                
                % set std
                obj.setVar( thisType, 'particulate', 'std', this_std);  
                
                % set binCount
                obj.setVar( thisType, 'particulate', 'binCount', binCount);
                
                % set all var/mean
                obj.setVar( thisType, 'particulate', 'all_var_over_mean', all_var_over_mean);
                obj.setVar( thisType, 'particulate', 'despiked_var_over_mean', despiked_var_over_mean);
                      
            end;   % end for loop a/c
            
            obj.L.info('ProcessedData.calcParticulate','End of Method');

         end; %#calcParticulate

        
        %% removeWLAfter750 & correctSpectralBandMismatch
        function obj = correctSpectralBandMismatch( obj, typeIn )
        %#removeWLAfter750 - Removes WL after 750 nm for c because the lookup 
        %#                 - table for scattering correction ends at 750
        %#                 - (a will be interpolated to c)
        %#                 - calls obj.cutWavelengths()
        %#correctSpectralBandMismatch - Deal with mismatch in spectral band positions between a and c measurements.
        %#                            - Interpolate a onto c, limiting range to where the two overlap
        %#
        %# SYNOPSIS correctSpectralBandMismatch( obj, typeIn )
        %# INPUT  obj    - the object
        %#        typeIn - 'a' or 'c'
        %# OUTPUT obj    - the object        
        %#
        %# SYNOPSIS removeWLAfter750( obj, typeIn )
        %# INPUT  obj    - the object
        %#        typeIn - "a" or "c"
        %# OUTPUT obj    - the object
        %#
             obj.L.info('ProcessedData.removeWLAfter750', 'Start of method');
              
             % parse input arguments
              if strcmp(typeIn, 'c')
                  obj.L.debug('ProcessedData.removeWLAfter750', 'c');
                  type = 'c';
              elseif strcmp(typeIn, 'a')
                 obj.L.debug('ProcessedData.removeWLAfter750', 'a');
                 type = 'a';
              else
                 obj.L.error('ProcessedData.removeWLAfter750', 'No appropriate type IN');
              end;
             
             % cp/ap
            thisType = sprintf('%sp', type);
            thisWL = sprintf('%sWavelengths', type);
             
              
            % get variables - if c passed in this is for cp
            wl = obj.meta.DeviceFile.(thisWL);
            data = obj.getVar('name', thisType,'data','data');
            uncertainty = obj.getVar('name', thisType, 'data', 'uncertainty');
            this_std = obj.getVar('name', thisType, 'data', 'std');
            binCount = obj.getVar('name', thisType, 'data', 'binCount');
            all_var_over_mean = obj.getVar('name', thisType, 'data', 'all_var_over_mean');
            despiked_var_over_mean = obj.getVar('name', thisType, 'data', 'despiked_var_over_mean');

            obj.L.debug('ProcessedData.removeWLAfter750', ...
                sprintf('size wl before: %u x %u', size(wl)));
            obj.L.debug('ProcessedData.removeWLAfter750', ...            
                sprintf('size data before: %u x %u', size(data)));
            obj.L.debug('ProcessedData.removeWLAfter750', ...
                sprintf('size uncertainty before: %u x %u', size(uncertainty)));
            obj.L.debug('ProcessedData.removeWLAfter750', ...
                sprintf('size std before: %u x %u', size(this_std)));
            
            % run procedure
            
            CUTOFF = 750;
            [rowIndex, ~] = find( wl > CUTOFF );
            
            % make a copy of the wavelengths
            if ~isempty(rowIndex)
                
                wlEndAt750 = wl;

                % for the new wavelengths, set the last wavelength over 750 to 750 and
                % blank the rest
                wlEndAt750(rowIndex(1)) = CUTOFF;
                wlEndAt750(rowIndex(2:end)) = [];
                
                obj.L.debug('ProcessedData.cutWavelengths', ...
                    sprintf('wlEndAt750 size: %u x %u', size(wlEndAt750)));
                
                % for the new 750 wavelength, find the appropriate values by interpolating
                % from the original data.  These all carry all 83 original wavelengths
                dataEndAt750 = interp1( wl, data', wlEndAt750, 'linear', 'extrap');
                dataEndAt750 = dataEndAt750';
                
                uncEndAt750 = interp1( wl, uncertainty', wlEndAt750, 'linear', 'extrap');
                uncEndAt750 = uncEndAt750';
                
                stdEndAt750 = interp1( wl, this_std', wlEndAt750, 'linear', 'extrap');
                stdEndAt750 = stdEndAt750';
                
                binCountEndAt750 = interp1( wl, binCount', wlEndAt750, 'linear', 'extrap');
                binCountEndAt750 = binCountEndAt750';
                
                all_var_over_meanEndAt750 = interp1( wl, all_var_over_mean', wlEndAt750, 'linear', 'extrap');
                all_var_over_meanEndAt750 = all_var_over_meanEndAt750';
                despiked_var_over_meanEndAt750 = interp1( wl, despiked_var_over_mean', wlEndAt750, 'linear', 'extrap');
                despiked_var_over_meanEndAt750 = despiked_var_over_meanEndAt750';
                
                if wl(1) > wlEndAt750(1)
                    obj.L.debug('ProcessedData.cutWavelengths',...
                        'First wavelength being interpolated AFTER first wl being interpolated to -- will be NaN');
                else
                    obj.L.debug('ProcessedData.cutWavelengths',...
                        'First wl being interpolated BEFORE or SAME as first WL being interpolated to');
                end;

             else   % if find > 750 returns nothing
                 
                 obj.L.debug('ProcessedData.cutWavelengths', 'find>750 returned nothing');
                 wlEndAt750 = wl;
                 dataEndAt750 = data;
                 uncEndAt750 = uncertainty;
                 stdEndAt750 = this_std;
                 binCountEndAt750 = binCount;
                 all_var_over_meanEndAt750 = all_var_over_mean;
                 despiked_var_over_meanEndAt750 = despiked_var_over_mean;
            end;
            
            % get variables we need:  c wavelengths, cp_uncorr
%             % want most recent cp - i.e. the ones just cut off @ 750
            cWavelengths = wlEndAt750;
            cData = dataEndAt750;
            cUncertainty = uncEndAt750;
            cSTD = stdEndAt750;
            cBinCount = binCountEndAt750;
            c_all_var_over_mean = all_var_over_meanEndAt750;
            c_despiked_var_over_mean = despiked_var_over_meanEndAt750;
            
            aWavelengths = obj.meta.DeviceFile.aWavelengths;
            aData = obj.getVar('name', 'ap', 'data', 'data');
            aUncertainty = obj.getVar('name', 'ap', 'data', 'uncertainty');
            aSTD = obj.getVar('name', 'ap', 'data', 'std');
            aBinCount = obj.getVar('name', 'ap', 'data', 'binCount');
            a_all_var_over_mean = obj.getVar('name', 'ap', 'data', 'all_var_over_mean');
            a_despiked_var_over_mean = obj.getVar('name', 'ap', 'data', 'despiked_var_over_mean');
            
            agData = obj.getVar('name', 'aFSW', 'data', 'interpolatedData');
            agUncertainty = obj.getVar('name', 'aFSW', 'data', 'interpolatedUncertainty');
            
            obj.L.debug('ProcessedData.correctSpectralBandMismatch', ...
                sprintf('size c wl before: %u x %u', size(cWavelengths)));
            obj.L.debug('ProcessedData.correctSpectralBandMismatch', ...            
                sprintf('size adata before: %u x %u', size(aData)));
            obj.L.debug('ProcessedData.correctSpectralBandMismatch', ...
                sprintf('size uncertainty before: %u x %u', size(aUncertainty)));
            obj.L.debug('ProcessedData.correctSpectralBandMismatch', ...
                sprintf('size awavelengths before: %u x %u', size(aWavelengths)));  
            obj.L.debug('ProcessedData.correctSpectralBandMismatch', ...
                sprintf('size aSTD before: %u x %u', size(aSTD)));  
            
            % run procedure (in this case interp1)
            
            % interpolate ap to cp wavelengths
            apMatchedData = interp1( aWavelengths, aData', cWavelengths, 'linear', 'extrap');
            apMatchedData = apMatchedData';
            apMatchedUnc = interp1( aWavelengths, aUncertainty', cWavelengths, 'linear', 'extrap');
            apMatchedUnc = apMatchedUnc';
            apMatchedSTD = interp1( aWavelengths, aSTD', cWavelengths, 'linear', 'extrap');
            apMatchedSTD = apMatchedSTD';
            apMatchedBinCount = interp1( aWavelengths, aBinCount', cWavelengths, 'linear', 'extrap');
            apMatchedBinCount = apMatchedBinCount';
            apMatchedAllVarOverMean = interp1( aWavelengths, a_all_var_over_mean', cWavelengths, 'linear', 'extrap');
            apMatchedAllVarOverMean  = apMatchedAllVarOverMean';
            apMatchedDespikedVarOverMean = interp1( aWavelengths, a_despiked_var_over_mean', cWavelengths, 'linear', 'extrap');
            apMatchedDespikedVarOverMean = apMatchedDespikedVarOverMean';
            
            agMatchedData = interp1( aWavelengths, agData', cWavelengths, 'linear', 'extrap');
            agMatchedData = agMatchedData';
            agMatchedUnc = interp1( aWavelengths, agUncertainty', cWavelengths, 'linear', 'extrap');
            agMatchedUnc = agMatchedUnc';
            
            % CALCULATING KEEP WAVELENGTHS INDEX, BUT NOT USING
            keepWavelengths = ~all(isnan(apMatchedData));
            obj.L.info('ProcessedData.correctSpectralBandMismatch', ...
                sprintf('keepWavlengths: %u x %u', size(keepWavelengths)));
                
            % SET VARS
            obj.setVar('ap', 'matchedWL', 'data', apMatchedData);
            obj.setVar('ap', 'matchedWL', 'uncertainty', apMatchedUnc);
            obj.setVar('ap', 'matchedWL', 'wavelengths', cWavelengths);
            obj.setVar('ap', 'matchedWL', 'std', apMatchedSTD);
            obj.setVar('ap', 'matchedWL', 'binCount', apMatchedBinCount);
            obj.setVar('ap', 'matchedWL', 'all_var_over_mean', apMatchedAllVarOverMean);
            obj.setVar('ap', 'matchedWL', 'despiked_var_over_mean', apMatchedDespikedVarOverMean);
            
            obj.setVar('aFSW', 'matchedWL', 'data', agMatchedData);
            obj.setVar('aFSW', 'matchedWL', 'uncertainty', agMatchedUnc);
            obj.setVar('aFSW', 'matchedWL', 'wavelengths', cWavelengths);
            
            obj.setVar('cp', 'matchedWL', 'data', cData);
            obj.setVar('cp', 'matchedWL', 'uncertainty', cUncertainty);
            obj.setVar('cp', 'matchedWL', 'wavelengths', cWavelengths);
            obj.setVar('cp', 'matchedWL', 'std', cSTD);
            obj.setVar('cp', 'matchedWL', 'binCount', cBinCount);
            obj.setVar('cp', 'matchedWL', 'all_var_over_mean', c_all_var_over_mean);
            obj.setVar('cp', 'matchedWL', 'despiked_var_over_mean', c_despiked_var_over_mean);
            
            % get c interpolated and rename as matched cg
            cgData = obj.getVar('name', 'cFSW', 'data', 'interpolatedData');
            cgUncertainty = obj.getVar('name', 'cFSW', 'data', 'interpolatedUncertainty');
            
            obj.setVar('cFSW', 'matchedWL', 'data', cgData);
            obj.setVar('cFSW', 'matchedWL', 'uncertainty', cgUncertainty);
            obj.setVar('cFSW', 'matchedWL', 'wavelengths', cWavelengths);
            
            obj.L.info('ProcessedData.correctSpectralBandMismatch', 'End of method');
           
        end  %#correctSpectralBandMismatch
         
        %% scatteringCorr
        function obj = scatteringCorr(obj, typeIn)
        %#scatteringCorr - performs residual temperature/scattering correction specified in input parameter 'typeIn'
        %#               - sets appropriate data
        %#               - 'Sullivan_etal_2006_instrumentspecific.xls'
        %#               - calls ResidTempScatCorr()
        %#               - calls getPsiT()
        %#
        %# SYNOPSIS myMethod(obj)
        %# INPUT obj    - the object
        %#       typeIn - the type of correction:  'SLADE', 'ROTTGERS', 'FLAT'
        %# OUTPUT obj   - the object
        %#
        

            obj.L.info('ProcessedData.scatteringCorr:', 'Start');
              
            % get vars
            % cWavelengths = obj.getVar('name', 'cp', 'data', 'wavelengths');
            apUncorr = obj.getVar('name', 'ap', 'data', 'data', 'level', 'matchedWL');
            cpUncorr = obj.getVar('name', 'cp', 'data', 'data', 'level', 'matchedWL');
            wavelengths = obj.getVar('name', 'cp', 'data', 'wavelengths', 'level', 'matchedWL');
              
            psiT = getPsiT('Sullivan_etal_2006_instrumentspecific.xls', wavelengths);
            psiT = psiT';

            
            % call function and set vars
            if strcmp(typeIn, 'SLADE')
              obj.L.debug('ProcessedData.scatteringCorr', 'SLADE');
              [ap_TSalScatCorr, ap_uncorr_ref, fiterr, deltaT] = ...
                  ResidTempScatCorr( apUncorr, cpUncorr, wavelengths, psiT, 'Slade');

              obj.setVar('ap', 'corrected', 'data_slade', ap_TSalScatCorr);
              obj.setVar('ap', 'corrected', 'uncorr_ref_slade', ap_uncorr_ref);
              obj.setVar('ap', 'corrected', 'wavelengths_slade', wavelengths);
              obj.setVar('ap', 'corrected', 'fiterr_slade', fiterr);
              obj.setVar('ap', 'corrected', 'deltaT_slade', deltaT);
              
             elseif strcmp(typeIn, 'ROTTGERS')
              obj.L.debug('ProcessedData.scatteringCorr', 'ROTTGERS');
              [ap_TSalScatCorr, ~, fiterr, deltaT] = ...
                  ResidTempScatCorr( apUncorr, cpUncorr, wavelengths, psiT, 'Rottgers');

              obj.setVar('ap', 'corrected', 'data_rottgers', ap_TSalScatCorr);
              obj.setVar('ap', 'corrected', 'wavelengths_rottgers', wavelengths);
              obj.setVar('ap', 'corrected', 'fiterr_rottgers', fiterr);
              obj.setVar('ap', 'corrected', 'deltaT_rottgers', deltaT);

            elseif strcmp(typeIn, 'FLAT')
              obj.L.debug('ProcessedData.scatteringCorr', 'FLAT');
              [ap_TSalScatCorr, ~, fiterr, deltaT] = ...
                  ResidTempScatCorr( apUncorr, cpUncorr, wavelengths, psiT, 'Flat');
              obj.setVar('ap', 'corrected', 'data_flat', ap_TSalScatCorr);
              obj.setVar('ap', 'corrected', 'wavelengths_flat', wavelengths);
              obj.setVar('ap', 'corrected', 'fiterr_flat', fiterr);
              obj.setVar('ap', 'corrected', 'deltaT_flat', deltaT);
            else
              obj.L.error('ProcessedData.scatteringCorr', 'No correct type IN');
            end;
            
            obj.L.info('ProcessedData.scatteringCorr:', 'End');
            
        end %#scatteringCorr
         
        %% attenuationCorr 
        function obj = attenuationCorr(obj, typeIn)
        %#attenuationCorr - calls external function AttTempCorr()
        %#                - calls getPsiT()
        %#
        %# SYNOPSIS attenuationCorr(obj, typeIn)
        %# INPUT  obj    - the object
        %#        typeIn - either "WITHA" or "WITHOUTA"
        %# OUTPUT obj    - the object
        %#             
            obj.L.info('ProcessedData.attenuationCorr:', 'Start');
             
            % get vars
            cpUncorr = obj.getVar('name', 'cp', 'data', 'data', 'level', 'matchedWL');
            wavelengths = obj.getVar('name', 'cp', 'data', 'wavelengths', 'level', 'matchedWL');
           
            psiT = getPsiT('Sullivan_etal_2006_instrumentspecific.xls', wavelengths);
            psiT = psiT';

            deltaTname = sprintf('deltaT_%s', lower(obj.meta.Params.PROCESS.SCATTERING_CORR));
            
             if strcmp(typeIn, 'WITHA')
                  obj.L.debug('ProcessedData.attenuationCorr', 'WITHA');
                  
                  % get deltaT from processing A
                  deltaT = obj.getVar('name', 'ap', 'data', deltaTname, 'level', 'corrected');
                  
                  [cpCorr, fiterr] = AttTempCorr(cpUncorr, wavelengths, psiT,'WITHA', deltaT);
                  obj.setVar('cp', 'corrected', 'data', cpCorr);
                  obj.setVar('cp', 'corrected', 'wavelengths', wavelengths);
                  obj.setVar('cp', 'corrected', 'fiterr', fiterr); 
                  
             elseif strcmp(typeIn, 'WITHOUTA')
                 obj.L.debug('ProcessedData.attenuationCorr', 'WITHOUTA');
                 
                 [cpCorr, fiterr] = AttTempCorr(cpUncorr, wavelengths, psiT,'WITHOUTA');
                 
                  obj.setVar('cp', 'corrected', 'data', cpCorr);
                  obj.setVar('cp', 'corrected', 'wavelengths', wavelengths);
                  obj.setVar('cp', 'corrected', 'fiterr', fiterr);                 
                 
                 
             else
                 obj.L.error('ProcessedData.attenuationCorr', 'No correct type IN');
             end;
             obj.L.info('ProcessedData.attenuationCorr:', 'End');
         end;

         
        %%  computeAPUncertaintyBetweenCorrections
        function obj = computeAPUncertaintyBetweenCorrections(obj)
        %#computeAPUncertaintyBetweenCorrections - calculate the uncertainty for ap
        %#                                       - by abs(slade -
        %rottgers/2) if both slade and rottgers calculations have been
        %done. 
        %#                                       - also calculates uncertainty
        %aTSWSTD./sqrt(obj.meta.Params.PROCESS.UNCERTAINTY_N); if mean used for
        %binning
        %#                                       - or calculates
        %uncertainty:
        %aTSWBinVariability./sqrt(obj.meta.Params.PROCESS.UNCERTAINTY_N); if median
        %used for binning
        %
        %#
        %# SYNOPSIS obj = computeAPUncertaintyBetweenCorrections(obj)
        %# INPUT obj: the object
        %# OUTPUT obj: the object
        %#
            obj.L.info('ProcessedData.computeAPUncertaintyBetweenCorrections:', 'Start');
             
            % get variables
            if obj.meta.Params.PROCESS.SCATTERING_CORR_SLADE && obj.meta.Params.PROCESS.SCATTERING_CORR_ROTTGERS

                slade = obj.getVar('name', 'ap', 'data', 'data_slade', 'level', 'corrected');
                rottgers = obj.getVar('name', 'ap', 'data', 'data_rottgers', 'level', 'corrected');
                
                %changed on 29-Jul-16
                uncertainty = abs(slade - rottgers)/2;
                
                % set variable:
                obj.setVar('ap', 'corrected', 'uncertainty_between_corrections', uncertainty);
            end;
            
            aInterpBinVariability = obj.getVar('name', 'aFSW', 'data', 'interpolatedUncertainty');

            % should already be calculated now:
            aBinUncertainty = obj.getVar('name', 'aTSW', 'data', 'uncertainty');
            uncertainty = aInterpBinVariability + aBinUncertainty;
            obj.setVar('ap', 'corrected', 'uncertainty', uncertainty);
          
             
            obj.L.info('ProcessedData.computeAPUncertaintyBetweenCorrections:', 'End');            
             
         end; %#computeAPUncertaintyBetweenCorrections

         
 %% unsmooth
        function unsmooth(obj, varargin)
        %#unsmooth - unsmooths a/c data
        %#         - calls spectral_unsmooth external function
        %#
        %# SYNOPSIS unsmooth(obj, varargin)
        %# INPUT  obj       - the object
        %#        dataType1 - the primary type of the data being binned, i.e. "a" or "c"
        %#                  - if no type is given, both will be run
        %# OUTPUT obj       - the object
        %# 
                         
            obj.L.info('ProcessedData.unsmooth()','Start of Method');

            unsmoothA = true;
            unsmoothC = true;

            % parse optional input parameters
            if (~isempty(varargin))
                
                % set to false initially
                unsmoothA = false;
                unsmoothC = false;
                
                for iArgs = 1:length(varargin)
                    switch varargin{iArgs}
                        case {'a'}
                            unsmoothA = true;
                        case {'c'}
                            unsmoothC = true;
                        otherwise
                            
                           obj.L.error('ProcessedData.unsmooth()', 'Invalid argument');
                    end   % switch
                end   % for
            else
                obj.L.debug('ProcessedData.unsmooth()', 'no args in');
            end;   % if varargin is empty
            
            
            
            if unsmoothA
                 % get ap - for each correction - more than one correction is possible
                if obj.meta.Params.PROCESS.SCATTERING_CORR_SLADE
                    ap_data_slade = obj.getVar('name', 'ap', 'data', 'data_slade', 'level', 'corrected');
                    wavelengths = obj.getVar('name', 'ap', 'data', 'wavelengths_slade', 'level', 'corrected');
                    unsmooth = spectral_unsmooth( wavelengths, ap_data_slade, 1);
                    obj.setVar('ap', 'unsmoothed', 'data_slade', unsmooth );
                end;
                if obj.meta.Params.PROCESS.SCATTERING_CORR_ROTTGERS
                    ap_data_rottgers = obj.getVar('name', 'ap', 'data', 'data_rottgers', 'level', 'corrected');
                    wavelengths = obj.getVar('name', 'ap', 'data', 'wavelengths_rottgers', 'level', 'corrected');
                    unsmooth = spectral_unsmooth( wavelengths, ap_data_rottgers, 1);
                    obj.setVar('ap', 'unsmoothed', 'data_rottgers', unsmooth );
                end;
                if obj.meta.Params.PROCESS.SCATTERING_CORR_FLAT
                    ap_data_flat = obj.getVar('name', 'ap', 'data', 'data_flat', 'level', 'corrected');
                    wavelengths = obj.getVar('name', 'ap', 'data', 'wavelengths_flat', 'level', 'corrected');
                    unsmooth = spectral_unsmooth( wavelengths, ap_data_flat, 1);
                    obj.setVar('ap', 'unsmoothed', 'data_flat', unsmooth );
                end;
             end %# if unsmoothA
             
            if unsmoothC

                % get cp
                cp_data = obj.getVar('name', 'cp', 'data', 'data', 'level', 'corrected');
                wavelengths = obj.getVar('name', 'cp', 'data', 'wavelengths', 'level', 'corrected');
                unsmooth = spectral_unsmooth( wavelengths, cp_data, 1);
                obj.setVar('cp', 'unsmoothed', 'data', unsmooth );  
                
            end %# if unsmoothC

            obj.L.info('ProcessedData.unsmooth','End of Method');
            
        end
        
        %% flagSuspectBins
        function obj = flagSuspectBins(obj, varargin)
        %calcSuspectData: sets variable 'Suspect Data' for FSW and TSW.  
        %Makes it easier to identify suspect data in plots, etc.
        %#
        %# SYNOPSIS obj = calcSuspectData(obj, varargin)
        %# INPUT obj: the object
        %#       dataType1      - the primary type of the data being binned, i.e. "a" or "c"
        %# OUTPUT obj: the object
        %#            
            obj.L.info('ProcessedData.flagSuspectBins','Start of Method');
            
            dataType1 = '';

            % parse optional input parameters
            if (~isempty(varargin))
                for iArgs = 1:length(varargin)
                    switch varargin{iArgs}
                        case {'a'}
                            dataType1 = {'a'};
                        case {'c'}
                            dataType1 = {'c'};
                        otherwise
                           obj.L.error('ProcessedData.flagSuspectBins', 'Invalid argument');
                    end   % switch
                end   % for
            else
                obj.L.debug('ProcessedData.flagSuspectBins', 'no args in');
                
                % by default, run both data types a & c
                dataType1 = {'a','c'};
                
            end;   % if varargin is empty
            
            % for a,c
            for iType1 = 1:length(dataType1)
                obj.L.info('ProcessedData.flagSuspectBins', sprintf('Starting %s', dataType1{iType1}));

                if strcmp(dataType1{iType1}, 'a')
                    stdThreshold = obj.meta.Params.PROCESS.a_STD_THRESH; %.015;
                else
                    stdThreshold = obj.meta.Params.PROCESS.c_STD_THRESH; %.03;
                end;
                
                % set the types of data that come before computing p
                thisType = sprintf('%sp', dataType1{iType1} );
                tswType = sprintf('%sTSW', dataType1{iType1} );
                fswType = sprintf('%sFSW', dataType1{iType1} );
 
                % 1.  get data                
                % get TSW bins and Interpolated FSW bins for calculating p
                % need: TSW.bin_median
                TSW_bin_median = obj.getVar('name', tswType, 'level', 'binned', 'data', 'median');
                % need: TSW.bin_mean
                TSW_bin_mean = obj.getVar('name', tswType, 'level', 'binned', 'data', 'mean');
                % need: FSW.interp_median
                FSW_interp_median = obj.getVar('name', fswType, 'level', 'filtered', 'data', 'interpolatedData');
                % need: TSW.bin_std
                TSW_bin_std = obj.getVar('name', tswType, 'level', 'binned', 'data', 'std');
                % need: rmoveIndex for timestamps beyond 24 hours
                removeIndex = obj.getVar('name', thisType, 'data', 'removeIndex');
            
                % 2.  set up flags
            
                binFlags = zeros(size(TSW_bin_median));
                binFlags(:,:) = 1; % 1 == good
            
                % 3.  Filter 1 - mean median filter
                mean_median_fail_index = ...
                    (abs(TSW_bin_median - TSW_bin_mean))./(TSW_bin_median-FSW_interp_median) ...
                    > max(0.3 , 0.001./(TSW_bin_median-FSW_interp_median));
                
                if sum(sum(mean_median_fail_index)) > 0                   
                    binFlags(mean_median_fail_index) = 3; % suspect

                    obj.L.info('ProcessedData.flagSuspectBins',...
                        sprintf('mean_median suspect bins: %u', ...
                        sum(sum(mean_median_fail_index))));
                    obj.L.info('ProcessedData.flagSuspectBins',...
                        sprintf('mean_median suspectspectra: %u', ...
                        sum(any(mean_median_fail_index,2))));                    
                end;

                % 4. Filter 2 - std threshold
                % if TSW.bin_std > max (stdThreshold, .5(TSW.bin_mean)
                % then 
                % flag as suspect
                std_fail_index = TSW_bin_std > stdThreshold;

                if sum(sum(std_fail_index)) > 0                   
                    binFlags(std_fail_index) = 3;  % 1 == good
                    obj.L.info('ProcessedData.flagSuspectBins',...
                        sprintf('std suspect bins: %u', ...
                        sum(sum(std_fail_index)))); 
                    obj.L.info('ProcessedData.flagSuspectBins',...
                        sprintf('std suspectspectra: %u', ...
                        sum(any(std_fail_index,2))));                       
                end;
            
                % 4.  remove extra rows from flags
                binFlags = removerows(binFlags, removeIndex);
                
                % 5. set flags
            
                obj.setVar(thisType, 'unsmoothed', 'flags', binFlags );  
                obj.L.info('ProcessedData.applyFlags',...
                    sprintf('total # bins: %u', ...
                    sum(any(binFlags == 3,2)))); 
             end;  %iType1 = 1:length(dataType1)
           

        obj.L.info('ProcessedData.flagSuspectBins','End of Method');                
        end  %#flagSuspectBins   
        

        function obj = removeSuspectBins(obj, varargin)
        %Makes it easier to identify suspect data in plots, etc.
        %#
        %# SYNOPSIS obj = calcSuspectData(obj, varargin)
        %# INPUT obj: the object
        %#       dataType1      - the primary type of the data being binned, i.e. "a" or "c"
        %# OUTPUT obj: the object
        %#            
            obj.L.info('ProcessedData.removeSuspectBins','Start of Method');

            % get data & flags
                
            ap_flags = obj.getVar('name','ap', 'level', 'unsmoothed', 'data', 'flags');
            cp_flags = obj.getVar('name','cp', 'level', 'unsmoothed', 'data', 'flags');
            
            % apply index to ther fields needed for printgin
             if obj.meta.Params.PROCESS.SCATTERING_CORR_SLADE
                obj.L.info('ProcessedData.removeSuspectBins','doing ap slade');

                ap_data_slade = obj.getVar('name', 'ap', 'data', 'data_slade'); 
                
                [good, suspect] = obj.applyFlags( ap_data_slade, ap_flags);
                obj.setVar( 'ap', 'flagsApplied', 'suspect_data_slade', suspect);
                obj.setVar( 'ap', 'flagsApplied', 'data_slade', good);
            end;
            if obj.meta.Params.PROCESS.SCATTERING_CORR_ROTTGERS
                ap_data_rottgers = obj.getVar('name', 'ap', 'data', 'data_rottgers'); 
                [good, suspect] = obj.applyFlags( ap_data_rottgers, ap_flags);
                obj.setVar( 'ap', 'flagsApplied', 'suspect_data_rottgers', suspect);
                obj.setVar( 'ap', 'flagsApplied', 'data_rottgers', good);
            end;
            if obj.meta.Params.PROCESS.SCATTERING_CORR_FLAT
                ap_data_flat = obj.getVar('name', 'ap', 'data', 'data_flat'); 
                [good, suspect] = obj.applyFlags( ap_data_flat, ap_flags);
                obj.setVar( 'ap', 'flagsApplied', 'suspect_data_flat', suspect);
                obj.setVar( 'ap', 'flagsApplied', 'data_flat', good);                
            end;
         
            cp_data = obj.getVar('name', 'cp', 'data', 'data'); %, 'level', 'corrected');
            [good, suspect] = obj.applyFlags( cp_data, cp_flags);
            obj.setVar( 'cp', 'flagsApplied', 'suspect_cp_data', suspect);
            obj.setVar( 'cp', 'flagsApplied', 'cp_data', good);  
            
            cp_std = obj.getVar('name', 'cp', 'data', 'std');
            [good, suspect] = obj.applyFlags( cp_std, cp_flags);
            obj.setVar( 'cp', 'flagsApplied', 'suspect_cp_std', suspect);
            obj.setVar( 'cp', 'flagsApplied', 'cp_std', good);            
            
            ap_std = obj.getVar('name', 'ap', 'data', 'std');
            [good, suspect] = obj.applyFlags( ap_std, ap_flags);
            obj.setVar( 'ap', 'flagsApplied', 'suspect_ap_std', suspect);
            obj.setVar( 'ap', 'flagsApplied', 'ap_std', good);
            
            cp_bin_count = obj.getVar('name', 'cp', 'data', 'binCount');
            [good, suspect] = obj.applyFlags( cp_bin_count, cp_flags);
            obj.setVar( 'cp', 'flagsApplied', 'suspect_cp_bin_count', suspect);
            obj.setVar( 'cp', 'flagsApplied', 'cp_bin_count', good);
            
            ap_bin_count = obj.getVar('name', 'ap', 'data', 'binCount');
            [good, suspect] = obj.applyFlags( ap_bin_count, ap_flags);
            obj.setVar( 'ap', 'flagsApplied', 'suspect_ap_bin_count', suspect);
            obj.setVar( 'ap', 'flagsApplied', 'ap_bin_count', good);            
                
        end;
        function [ goodDataOut, suspectDataOut ] = applyFlags( obj, dataIn, flagsIn )
            obj.L.info('ProcessedData.applyFlags','Start of Method');
            obj.L.info('ProcessedData.applyFlags',...
                    sprintf('total # flags: %u', ...
                    sum(any(flagsIn == 3,2)))); 
            flags = flagsIn;
            all_data = dataIn;
            suspectBinMatrixIndex(:,:) = (flags(:,:) == 3);
            % set to row index, not matrix index
            suspectBinIndex(:,:) = any(suspectBinMatrixIndex, 2);
            numBadBins = sum(suspectBinIndex);
            obj.L.info('ProcessedData.applyFlags',...
                sprintf('num bad bins being blanked: %u', numBadBins));
            suspectBinMatrixIndex(suspectBinIndex,:) = 1;
            suspect_data = zeros(size(all_data));
            suspect_data(:) = NaN;
            suspect_data(suspectBinMatrixIndex(:,:)) = ...
                all_data(suspectBinMatrixIndex(:,:));
            
            good_data = zeros(size(all_data));
            good_data(:) = NaN;
            good_data(~suspectBinMatrixIndex(:,:)) = all_data(~suspectBinMatrixIndex(:,:));
            
            goodDataOut = good_data;
            suspectDataOut = suspect_data;
            obj.L.info('ProcessedData.applyFlags','End of Method');

        end;

%%
    %-----------------------------------------------------------------
    % PLOTS
    %-----------------------------------------------------------------
        
        
        function plotSuspectData(obj, fignum, wavelengthToPlot)
            % plot original data against "good" and "suspect" FSW data, marking bins
            figure(fignum);
            
            ax1 = subplot(2,1,1);
            hold on;
            grid on;
            plot(obj.var.a.L1.timestamps, obj.var.a.L1.data(:,wavelengthToPlot), 'c');
            scatter(obj.var.ap.L5.timestamps, obj.var.ap.L9.data_slade(:, wavelengthToPlot), 'ko');
            scatter(obj.var.ap.L5.timestamps, ...
                obj.var.ap.L9.suspect_data_slade(:, wavelengthToPlot), 'ro');
            % bug fix for matlab legend not showing different colors in
            % scatter:
            [h,~] = legend('raw data','good ap bins', 'suspect ap bins');           
            title('ap data - good vs. suspect bins')
            dynamicDateTicks;

            ax2 = subplot(2,1,2);
            hold on;
            grid on;
            plot(obj.var.c.L1.timestamps, obj.var.c.L1.data(:, wavelengthToPlot), 'c');
            scatter(obj.var.cp.L5.timestamps, obj.var.cp.L9.cp_data(:, wavelengthToPlot), 'ko');
            scatter(obj.var.cp.L5.timestamps, ...
                obj.var.cp.L9.suspect_cp_data(:, wavelengthToPlot), 'ro');
            title('cp data - good vs. suspect bins');
            % bug fix for matlab legend not showing different colors in
            % scatter:            
            [h,~] = legend('raw data', 'good cp bins', 'suspect cp bins');
            dynamicDateTicks;
            linkaxes([ax1, ax2], 'x')
        end

        function plotACInterpolatedData(obj, fignum, wavelengthToPlot)
            figure(fignum);
            ax1 = subplot(2,1,1);
            hold on;
            grid on;
            plot(obj.var.c.L1.timestamps, obj.var.c.L1.data(:, wavelengthToPlot), 'y');
            scatter(obj.var.cFSW.L3.binnedTime, obj.var.cFSW.L4.interpolatedSectionMedians(:,20),'o','MarkerEdgeColor','blue');
            scatter(obj.var.cFSW.L4.binnedTime, obj.var.cFSW.L4.interpolatedData(:,20),10,'o','fill','MarkerEdgeColor','green');
            legend('Synchronized c data', 'FSW Bins', 'Interpolated data');
            title('Interpolated Filtered c Data - using median','fontsize',12);

            dynamicDateTicks;
            ax2 = subplot(2,1,2);
            hold on;
            grid on;
            plot(obj.var.a.L1.timestamps, obj.var.a.L1.data(:, wavelengthToPlot), 'y');
            scatter(obj.var.aFSW.L3.binnedTime, obj.var.aFSW.L4.interpolatedSectionMedians(:,20),'o','MarkerEdgeColor','blue');
            scatter(obj.var.aFSW.L4.binnedTime, obj.var.aFSW.L4.interpolatedData(:,20),10,'o','fill','MarkerEdgeColor','green');
            title('Interpolated Filtered a Data - using median','fontsize',12);
            legend('Synchronized a data', 'FSW Bins', 'Interpolated data');
            dynamicDateTicks;
%             linkaxes([ax1, ax2], 'xy')
        end       
         
        function plotCpVsTime(obj, fignum)
            figure(fignum)
            ax1 = subplot(3,1,1);
            hold on;
            grid on;
            plot(obj.var.cFSW.L4.binnedTime, obj.var.cFSW.L4.interpolatedData)
            legend('interpolated filtered')
            dynamicDateTicks;
            
            ax2 = subplot(3,1,2);
            hold on;
            grid on;
            plot(obj.var.cTSW.L3.binnedTime, obj.var.cTSW.L3.median)
            legend('cTSW bins')
            dynamicDateTicks;
            
            ax3 = subplot(3,1,3);
            hold on;
            grid on;
            plot(obj.var.cp.L5.timestamps, obj.var.cp.L5.data)
            legend('cp')
            dynamicDateTicks;
            
            linkaxes([ax1,ax2,ax3], 'x')
        end
        
        function plotApVsTime(obj, fignum)
            
            figure(fignum)
            ax1 = subplot(3,1,1);
            hold on;
            grid on;
            plot(obj.var.aFSW.L4.binnedTime, obj.var.aFSW.L4.interpolatedData)
            legend('interpolated filtered')
            dynamicDateTicks;
            
            ax2 = subplot(3,1,2);
            hold on;
            grid on;
            plot(obj.var.aTSW.L3.binnedTime, obj.var.aTSW.L3.median)
            legend('aTSW bins')
            dynamicDateTicks;
            
            ax3 = subplot(3,1,3);
            hold on;
            grid on;
            plot(obj.var.ap.L5.timestamps, obj.var.ap.L5.data)
            legend('ap')
            dynamicDateTicks;
            
            linkaxes([ax1,ax2,ax3], 'x')
            
        end
        
%         function plotCpVsWL(obj, fignum, level)
%             figure(fignum)
%             plot(obj.meta.DeviceFile.cWavelengths, obj.var.cp.(level).uncorr(spectraToPlot,:));
%         end
%         
%         function plotApVsWL(obj, fignum, level)
%             figure(fignum)
%             plot(obj.var.ap.(level).wavelengths, obj.var.ap.(level).data);
%         end
        
        function plotApCorr(obj, fignum, correctionMethod)
            % correctionMethod = '
            figure(fignum)
            wavelengths = sprintf('wavelengths_%s', correctionMethod);
            data = sprintf('data_%s', correctionMethod);
            plot( obj.var.ap.L7.(wavelengths), obj.var.ap.L7.(data));
            legend(sprintf('ap - corrected - %s', correctionMethod));
        end
        function plotCpCorr( obj, fignum)
            figure(fignum)
            hold on;
            grid on;
            plot( obj.var.cp.L7.wavelengths, obj.var.cp.L7.data);
            legend('cp - corrected - %s');
        end;

                
        function plotApWithoutSuspectData(obj, fignum, correctionMethod)
            % correctionMethod = '
            figure(fignum)
            wavelengths = sprintf('wavelengths_%s', correctionMethod);
            data = sprintf('data_%s', correctionMethod);
            plot( obj.var.ap.L7.(wavelengths), obj.var.ap.L9.(data));
            legend(sprintf('ap - corrected - %s', correctionMethod));
        end
        function plotCpWithoutSuspectData( obj, fignum)
            figure(fignum)
            hold on;
            grid on;
            plot( obj.var.cp.L7.wavelengths, obj.var.cp.L9.cp_data);
            legend('cp - corrected ');
        end;
            
        function plotCpCorrVsUncorr( obj, fignum )
            figure(fignum)
            tsToUse = find( ~isnan(obj.var.cp.L7.data(:,1)), 1, 'first');
            hold on; 
            grid on;
            plot( obj.var.cp.L7.wavelengths, obj.var.cp.L7.data(tsToUse,:), '*b')
            plot( obj.var.cp.L6.wavelengths, obj.var.cp.L6.data(tsToUse,:), '*c')
            legend('corrected cp', 'uncorrected cp')
            title('correcting c using a')
        end;
         
        function plotUnsmoothVsSmooth(obj, fignum, spectraToPlot, cSpectraToPlot )
            figure(fignum)
            ax1 = subplot(2,1,1);
            grid on;
            plot(obj.var.ap.L7.wavelengths_slade, obj.var.ap.L8.data_slade(spectraToPlot,:), 'b')
            hold on;
            plot(obj.var.ap.L7.wavelengths_slade, obj.var.ap.L7.data_slade(spectraToPlot,:), 'k')
            legend('ap - after unsmoothed', 'ap - before unsmoothing');
            
            ax2 = subplot(2,1,2);
            hold on;
            grid on;
            plot(obj.var.cp.L7.wavelengths, obj.var.cp.L8.data(cSpectraToPlot,:), 'b')
            plot(obj.meta.DeviceFile.cWavelengths, obj.var.cp.L5.data(cSpectraToPlot,:), 'k'); 
            legend('cp - after unsmoothed', 'cp - before unsmoothing');
        end;        
    % --------------------------------------------------------------------
    % END METHODS
    % --------------------------------------------------------------------
    

    
    end   %# methods
end   %# classdef
            
            
            