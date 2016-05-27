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
   end
   properties (Access = private)
       L   % logger
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
        
      
        
        function plotData(obj)
            
            obj.L.debug('AncillaryData.plotData', 'Start');
            ts = obj.DataObject;
            scatter(obj.DataObject.Time, obj.DataObject.Data, '.')
            dynamicDateTicks();
            
        end   % end plotData
        
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
            
            
            obj.L.debug('ProcessedData.setVar()', 'in method');

            level = obj.levelsMap(levelNameIn);
            obj.var.(varNameIn).(level).(dataNameIn) = dataIn;
            
            %find index number for this level
            levelInt = level(:,2);  % get just the number off
            levelInt = str2num(levelInt);
           

            if isempty(obj.levelsFlags)
                
                % if no flags exist at all yet -- create flags and index
                obj.L.debug('ProcessedData.setVar()', 'have no flags at all');
                levelFlags = 1:9;
                levelIndex = zeros(size(levelFlags));
                levelIndex = logical(levelIndex);
                obj.levelsFlags.(varNameIn).levelFlags = levelFlags;
                obj.levelsFlags.(varNameIn).levelIndex = levelIndex;
                obj.levelsFlags.(varNameIn).levelIndex(levelInt) = true;
                
            else
                obj.L.debug('ProcessedData.setVar()','we do have flags');
                
                if ~isfield(obj.levelsFlags, varNameIn )
                    
                    % it's not empty - but need to check if we have this level
                    obj.L.debug('ProcessedData.setVar()','have flags but not for this data');
                    levelFlags = 1:9;
                    levelIndex = zeros(size(levelFlags));
                    levelIndex = logical(levelIndex);
                    obj.levelsFlags.(varNameIn).levelFlags = levelFlags;
                    obj.levelsFlags.(varNameIn).levelIndex = levelIndex;
                    obj.levelsFlags.(varNameIn).levelIndex(levelInt) = true;
                    
                else
                    
                    % just update correct level
                    obj.L.debug('ProcessedData.setVar()','have flags for this data');
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
            
            obj.L.debug('ProcessedData.getVar()', 'in method');
            
            % check inputs
            % check varName is one that exists. Create if it doesn't exist?
            varNameIn = '';
            levelNameIn = '';
            dataNameIn = '';
            
            % loop through name/value pairs
            if (~isempty(varargin))
                iArg = 1;
                while iArg < nargin
%                     varargin{iArg};
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
                    obj.L.debug('ProcessedData.getVar()', ...
                        sprintf('level: %s', level));
                else
                    obj.L.error('ProcessedData.getVar', 'invalid data name');
                end;
                
            elseif ~isempty(varNameIn) && ~isempty(levelNameIn) && isempty(dataNameIn)
                obj.L.debug('ProcessedData.getVar()', ...
                    'have both name and level, don''t need data -- getting specific level');
                
                % lookup level
                level = obj.levelsMap(levelNameIn);

                if isfield(obj.var.(varNameIn), level )
                    varOut = obj.var.(varNameIn).(level);
                else
                    obj.L.error('ProcessedData.getVar', 'invalid level');
                end;

            elseif ~isempty(varNameIn) && isempty(levelNameIn)  && isempty(dataNameIn)
                obj.L.debug('ProcessedData.getVar()', ...
                    'have name, don''t need data -- need to find most recent level');
                
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
                                obj.L.debug('ProcessedData.getVar()', ...
                                    sprintf('setting data to level: %s', level));

                            else
                                obj.L.debug('ProcessedData.getVar()', ...
                                    'data doesn''t exist at this level');
                            end;
                        else
                            obj.L.debug('ProcessedData.getVar()', ...
                                'level doesn''t exist');
                        end;
                    end;   % for loop through levels
                else
                    obj.L.error('ProcessedData.getVar()', 'no data for this variable');
                end;  %isempty(thisLevels)              
            else
                obj.L.error('ProcessedData.getVar', 'problem');
            end; %~isempty(varNameIn) && ~isempty(levelNameIn) && ~isempty(dataNameIn)
                
        end; %#getVar
   end
     
end