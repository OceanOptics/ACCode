%%ACData Class
% 5/2016 - changed the way Transitions are handled to use a structure.
%        - change to inherit from AncillaryInlineData, not just
%        AncillaryData

classdef ACData < AncillaryInlineData
    properties
        Name
        Type = 'ac';
        Units
        
        % a timeseries object
        DataObject
        SyncDataObject
        
        % preprocessing variables
        SmoothData
        TSWIndex;
        FSWIndex;

        runningFSWmedian
        runningTSWmedian
        
        levelsMap         %# a map of the level names and numbers
        levelsFlags       %# a set of flags showing which levels are in use
        L
        var;        
    end
    properties (Access = private)
        Debug = 0;
    end
    methods
        
        function obj = ACData(nameIn, dataValuesIn, timestampsIn, unitsIn)
            %%% Pre-initialization %%%
            % Any code not using output argument (obj)
            
            %%% Object Initialization %%%
            % Call superclass constructor before accessing object
            % This statment cannot be conditionalized
                
            obj = obj@AncillaryInlineData(nameIn, dataValuesIn, timestampsIn, unitsIn);
            
            %%% Post-initialization %%%
            % Any code, including access to the object
            % adding variables to timeseries object.
            obj.L = log4m.getLogger();
            obj.L.debug('ACData.ACData()','Created object');
        end
        
        %% ----------------------------------------------------------------
        % setSmoothData
        % -----------------------------------------------------------------
        function setSmoothData(obj)
            obj.L.info('ACData.setSmoothData', 'start method');
            obj.SmoothData = filtfilt( ones(1,100)/100, 1, obj.DataObject.Data);
            obj.L.debug('ACData.setSmoothData', ...
                sprintf('size of orig data: %u x %u',size(obj.DataObject.Data)));
            obj.L.debug('ACData.setSmoothData', ...
                sprintf('size of  smooth data: %u x %u',size(obj.SmoothData)));
            obj.L.info('ACData.setSmoothData', 'end method');
        end

        %% ----------------------------------------------------------------
        % setRunningMedians
        % calculate the running medians - used to identify TSW and FSW
        % -----------------------------------------------------------------
        function setRunningMeds(obj, PPTimespan, SamplingFreq)
            obj.L.info('ACData.setRunningMeds', 'start method');
            [obj.runningFSWmedian, obj.runningTSWmedian] = ...
                findFSWTSWRunMeds(obj.SmoothData(:,20), obj.DataObject.Time, PPTimespan, SamplingFreq);
            obj.L.info('ACData.setRunningMeds', 'end method');
        end
         
        %% ----------------------------------------------------------------
        % findTransitions - NEW
        % -----------------------------------------------------------------
        function findAndSetTransitionPoints( obj, typeIn )
            obj.L.info('ACData.findAndSetTransitionPoints', 'start method');
            type = typeIn;
            timestamps = obj.DataObject.Time;
            
            % get opposite type running median to use
            if strcmp(type, 'TSW')
                runningMedian1 = obj.runningFSWmedian;
                runningMedian2 = obj.runningTSWmedian;
                
            elseif strcmp(type, 'FSW')
                runningMedian1 = obj.runningTSWmedian;
                runningMedian2 = obj.runningFSWmedian;
            else
                obj.L.error('ACData.findAndSetTransitionPoints', 'ERROR - incorrect type');
            end;
            transitions.(type) = findTransitionPoints( obj, runningMedian1, ...
                runningMedian2, timestamps);
            name = strcat(type, 'transitions');

            obj.setVar( 'TSW', 'raw', name, transitions.(type) );

            obj.L.info('ACData.findAndSetTransitionPoints', 'end method');
        end % function findAndSetTransitionPoints
        
        function transitions = findTransitionPoints( obj, runningMedian1In, ...
                runningMedian2In, timestampsIn )
            obj.L.info('ACData.findTransitionPoints', 'start method');
            
            runningMedian1 = runningMedian1In;
            runningMedian2 = runningMedian2In;
            timestamps = timestampsIn;
            
            nanArray = isnan(runningMedian1);
            rm = smooth(nanArray,151);

            % get rid of rounding error
            rm_round = (floor(rm*100))/100;
            rm = rm_round;

            % create index of values where it is just 1/filtered, not transition
            idxFiltered = (rm == 1);

            diffArray =  idxFiltered(2:end) - idxFiltered(1:end-1);
            diffArray(end+1) = NaN;

            startsIndex = find(diffArray > 0);
            endsIndex = find(diffArray < 0);
            for i = 1:length(startsIndex)
                obj.L.debug('ACData.findTransitionPoints',...
                    sprintf('startsIndex(i)+76: %u', startsIndex(i) + 76));
                if startsIndex(i) + 76 > length(runningMedian2)
                    obj.L.error('ACData.findTransitionPoints',...
                        sprintf('About to run over beginning of runningmedian; i: %u', i));
                    % make this the end
                    startData(i) = runningMedian2(startsIndex(1));
                    startTime(i) = timestamps(startsIndex(1));                   
                else
                    obj.L.debug('ACData.findTransitionPoints','Starts index + 76 within end');
                    startData(i) = runningMedian2(startsIndex(i) + 76);
                    startTime(i) = timestamps(startsIndex(i) + 76);

                end
            end  % end startsIndex loop
            for i = 1:length(endsIndex)
                obj.L.debug('ACData.findTransitionPoints',...
                    sprintf('endsIndex(i)+76: %u', endsIndex(i) + 76));
                if endsIndex(i) - 76 < 1
                    obj.L.error('ACData.findTransitionPoints',...
                        sprintf('About to run over end of runningmedian; i: %u', i));
                    % make this beginning 
                     endData(i) = runningMedian2(endsIndex(end));
                     endTime(i) = timestamps(endsIndex(end));
                else
                     obj.L.debug('ACData.findTransitionPoints','ends index - 76 within beginning'); 
                     endData(i) = runningMedian2(endsIndex(i) - 76);
                     endTime(i) = timestamps(endsIndex(i) - 76);
                end
            end  % end endsIndex loop

            transitions.StartTime = startTime';
            transitions.EndTime = endTime';
            transitions.StartData = startData';
            transitions.EndData = endData';
            
            obj.L.info('ACData.findTransitionPoints', 'start method');            
            
        end;   % end function findTransitionPoints();
          
            
        %% ----------------------------------------------------------------
        % getTransitions(obj)
        %# SYNOPSIS obj = processBins(obj, varargin)
        %# INPUT  obj            - the object
        %#        dataType1      - the primary type of the data being binned, i.e. "a" or "c"
        %#        dataType2      - the secondary type of the data being binned, i.e. TSW or FSW
        %# OUTPUT obj            - the object
        %# 
        function [Transitions] = getTransitions(obj, levelIn, varargin)
            
            obj.L.info('getInitTransitions','Start of Method');
           
            % check varargin
            dataType1 = '';
            dataType2 = '';
            
            transitionLevel = levelIn;

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
                            
                           obj.L.error('getInitTransitions', 'Invalid argument');
                    end   % switch
                end   % for
            else
                obj.L.debug('getInitTransitions', 'no args in');
                dataType1 = {'a','c'};
                dataType2 = {'TSW','FSW'};
            end;   % if varargin is empty
            
            for iType1 = 1:length(dataType1)
                  for iType2 = 1:length(dataType2)

                      obj.L.debug('getInitTransitions', sprintf('type: %s', dataType1{iType1}));
                      obj.L.debug('getInitTransitions', sprintf('type: %s', dataType2{iType2}));
                      
                      thisType = sprintf('%s%s', dataType1{iType1}, dataType2{iType2});
                      transType = sprintf('%stransitions', dataType2{iType2});
                      thisTransitions = obj.getVar('name', 'TSW', ...
                          'level', transitionLevel, 'data', transType);
                      Transitions.(transType) = thisTransitions;
                  end;
            end;
            obj.L.info('getInitTransitions','End of Method');
        end;


        %% ----------------------------------------------------------------
        % setGoodTransitions
        % -----------------------------------------------------------------
        function setGoodTransitions(obj, goodStarts, goodEnds, type, level)
            
            name = strcat(type, 'transitions');
            goodStarts = goodStarts(~isnan(goodStarts));
            goodEnds = goodEnds(~isnan(goodEnds));
               
            [~,ind1] = min(abs(bsxfun(@minus,datenum(obj.DataObject.Time),datenum(goodStarts)')));
            % subset the data to the good times
            goodStartData = obj.DataObject.Data(ind1,20); 

            [~,ind1] = min(abs(bsxfun(@minus,datenum(obj.DataObject.Time),datenum(goodEnds)')));
            goodEndData = obj.DataObject.Data(ind1,20);

            transitions.StartTime = goodStarts;
            transitions.EndTime = goodEnds;
            transitions.StartData = goodStartData;
            transitions.EndData = goodEndData;

            obj.setVar( 'TSW', level, name, transitions );
       
        end
       
        %% ---------------------------------------------------------------
        % ----------------------------------------------------------------
        % Sync to should be run AFTER checkTransitions()
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        function syncTo(obj, AncillaryDataIn)
            obj.L.debug('ACData.syncTo', 'Start of Method');
            
            % WE ARE MATCHING FILTERED STARTS HERE, NOT TOTAL
            Transitions = obj.getTransitions('preprocessed');
            goodTrans = Transitions.FSWTransitions.StartTime;
            
            obj.L.debug('ACData.syncTo', ...
                sprintf('first good transition in FSW: %s', datestr(goodTrans(1))));
            
            if length(goodTrans) >= 1 
                obj.L.debug('ACData.syncTo', 'at least one goodTrans');
                dataInGood = AncillaryDataIn.TransStartTime(AncillaryDataIn.TransStartFlag == 1);
                [~,ind1] = min(abs(datenum(dataInGood) - datenum(goodTrans(1))));
                closestGoodTime = dataInGood(ind1,:);
                obj.L.debug('ACData.syncTo', ...
                    sprintf('closestGoodTime from ancillary: %s', datestr(closestGoodTime)));
                
                offset = closestGoodTime-goodTrans(1);
                offsetAsDateVec = datevec(offset);
                size(offsetAsDateVec)
                obj.L.debug('ACData.syncTo', ...
                    sprintf('offset: H: %u M: %u S: %u', offsetAsDateVec(4), offsetAsDateVec(5), offsetAsDateVec(6)));
                
                
                obj.L.debug('ACData.syncTo',...
                    sprintf('original data object start time: %s, end time: %s', ...
                    datestr(obj.DataObject.Time(1)), datestr(obj.DataObject.Time(end))));
                
                obj.SyncDataObject = obj.DataObject; 
                obj.SyncDataObject.Time = obj.SyncDataObject.Time + offset;
                
                obj.L.debug('ACData.syncTo',...
                    sprintf('synced data object start time: %s, end time: %s', ...
                    datestr(obj.SyncDataObject.Time(1)), datestr(obj.SyncDataObject.Time(end))));
                
               
                Transitions.TSWtranstions.StartTime = ...
                    Transitions.TSWtransitions.StartTime + offset;
                Transitions.TSWtransitions.EndTime = ...
                    Transitions.TSWtransitions.EndTime + offset;                
                Transitions.FSWTransitions.StartTime = ...
                    Transitions.FSWTransitions.StartTime + offset;
                Transitions.FSWTransitions.EndTime = ...
                    Transitions.FSWTransitions.EndTime + offset;  
                
                obj.setVar( 'TSW', 'preprocessed', 'TSWtransitions', ...
                    Transitions.TSWtransitions);
                obj.setVar( 'TSW', 'preprocessed', 'FSWtransitions', ...
                    Transitions.FSWTransitions);
               
            else
                obj.L.error('ACData.syncTo', 'No good transitions to sync to!');
            end
            
            obj.L.debug('ACData.syncTo', 'Start of Method');
        end
               
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        % offsetFSWStartTimes() 
        % This function will shift forward the FSW start times by a given
        % time (in seconds)
        function offsetFSWStartTimes(obj, offsetTimeIn)
        
            % NEED CHECK HERE FOR SYNCDATAOBJECT
            
            
            
            Transitions = obj.getTransitions('preprocessed');
            startTimes = Transitions.FSWtransitions.StartTime;
            
            % convert time in seconds to a date vec
            offsetDatenum = datenum(0, 0, 0, 0, 0, offsetTimeIn);
            % for each start time:
            
            for iStartTime = 1:numel(startTimes)
            %   add offsetTime 
                obj.L.debug('ACData.offsetFSWStartTimes', ...
                    sprintf('old time: %s', datestr(startTimes(iStartTime))));
                
                tempStartTime = startTimes(iStartTime) + offsetDatenum;
                
                obj.L.debug('ACData.offsetFSWStartTimes', ...
                    sprintf('temp start time: %s', datestr(tempStartTime)));  
                
                %   find index of new start time
                [~,ind1] = min(abs(datenum(obj.SyncDataObject.Time) - datenum( tempStartTime )));

                %   reasssign startData with correct time and data
                Transitions.FSWtransitions.StartTime(iStartTime) = ...
                    obj.SyncDataObject.Time(ind1,:);
                
                % just use col 20 here -- that's what orig running med used
                Transitions.FSWtransitions.StartData(iStartTime) = ...
                    obj.SyncDataObject.Data(ind1,20);                

            end   % end for loop
          obj.setVar( 'TSW', 'preprocessed', 'TSWtransitions', ...
                    Transitions.TSWtransitions);
            obj.setVar( 'TSW', 'preprocessed', 'FSWtransitions', ...
                    Transitions.FSWtransitions);      
        end  % end function offsetFSWStartTimes
        
        %% ----------------------------------------------------------------
        % -----------------------------------------------------------------
        % create two indexes TSWIndex, FSWIndex to mark TSW and FSW and
        % then flag transition periods as 'suspect' (3)
        function separateTSWFSW(obj)
            % create two NaN arrays the size of data
            
            % CHANGE THESE DATAOBJECTS BACKINTO DATASYNCOBJECTS
            % OR SUB DATASYNC FOR DATA
            
            Transitions = obj.getTransitions('preprocessed');
            FSWStartTimes = Transitions.FSWtransitions.StartTime;
            FSWEndTimes = Transitions.FSWtransitions.EndTime;
            TSWStartTimes = Transitions.TSWtransitions.StartTime;
            TSWEndTimes = Transitions.TSWtransitions.EndTime;
            
            %obj.TSWIndex = zeros( size( obj.SyncDataObject.Time ));
            % want 0s for false, not nans
            obj.FSWIndex = zeros( length( obj.SyncDataObject.Time ), 1);
            obj.TSWIndex = zeros( length( obj.SyncDataObject.Time ), 1);

            % get index of start times in data
            [~,indFSWStarts] = min(abs(bsxfun(@minus,...
                datenum(obj.SyncDataObject.Time),...
                datenum(FSWStartTimes)')));
            [~,indFSWEnds] = min(abs(bsxfun(@minus,...
                datenum(obj.SyncDataObject.Time),...
                datenum(FSWEndTimes)')));
             
            obj.L.debug('ACData.separateTSWFSW', ...
                sprintf('length of indFSWStarts: %u', length(indFSWStarts)));
            obj.L.debug('ACData.separateTSWFSW', ...
                sprintf('length of indFSWEnds: %u', length(indFSWEnds)));

            for iFSW = 1:length(indFSWStarts)
                obj.L.debug('ACData.separateTSWFSW', ...
                sprintf('marking between %u and %u with 1', indFSWStarts(iFSW), indFSWEnds(iFSW)));     
            
                obj.FSWIndex(indFSWStarts(iFSW):indFSWEnds(iFSW)) = 1;
            end;
            
            % convert index to logical type
            obj.FSWIndex = logical(obj.FSWIndex);

            obj.L.debug('ACData.separateTSWFSW', ...
                sprintf('FSW index min: %s', ...
                datestr(min(obj.SyncDataObject.Time(obj.FSWIndex == 1)))));
            obj.L.debug('ACData.separateTSWFSW', ...
                sprintf('FSW index max: %s', ...
                datestr(max(obj.SyncDataObject.Time(obj.FSWIndex == 1)))));  
            
            % ------------------------------------------------------------
            % Do same for TSW:
            % get index of start times in data
            [~,indTSWStarts] = min(abs(bsxfun(@minus,...
                datenum(obj.SyncDataObject.Time),...
                datenum(TSWStartTimes)')));
            [~,indTSWEnds] = min(abs(bsxfun(@minus,...
                datenum(obj.SyncDataObject.Time),...
                datenum(TSWEndTimes)')));
             
            obj.L.debug('ACData.separateTSWFSW', ...
                sprintf('length of indTSWStarts: %u', length(indTSWStarts)));
            obj.L.debug('ACData.separateTSWFSW', ...
                sprintf('length of indTSWEnds: %u', length(indTSWEnds)));

            for iTSW = 1:length(indTSWStarts)
             obj.L.debug('ACData.separateTSWFSW', ...
                sprintf('marking between %u and %u with 1', indTSWStarts(iTSW), indTSWEnds(iTSW)));                   
                obj.TSWIndex(indTSWStarts(iTSW):indTSWEnds(iTSW)) = 1;
            end;
            
            % convert index to logical type
            obj.TSWIndex = logical(obj.TSWIndex);

            obj.L.debug('ACData.separateTSWFSW', ...
                sprintf('TSW index min: %s', ...
                datestr(min(obj.SyncDataObject.Time(obj.TSWIndex == 1)))));
            obj.L.debug('ACData.separateTSWFSW', ...
                sprintf('TSW index max: %s', ...
                datestr(max(obj.SyncDataObject.Time(obj.TSWIndex == 1)))));  
           
        end   % end function
        
        %%
        % ----------------------------------------------------------------
        % plots
        % ----------------------------------------------------------------
        function plotData(obj)
            if obj.Debug == 1;
                disp('Stub for plotData() in ACData:');
            end;
            
            ts = obj.DataObject;
            plot(ts.Time, ts.Data(:,20), ':bo', ...
                'MarkerSize', 3, ...
                'MarkerFaceColor', 'b');
            
            xlabel('Timestamp')
            ylabel(obj.Name)
            title(obj.Name)
            dynamicDateTicks
        end
        function plotSmoothData(obj)
            if obj.Debug == 1;
                disp('plotSmoothData() in ACData');
            end;
            plot(obj.DataObject.Time, obj.SmoothData(:,20), 'go', ...
                'LineWidth', 1, ...
                'Color', 'b',...
                'MarkerSize', 3, ...
                'MarkerFaceColor', 'k')
            xlabel('Timestamp')
            ylabel(strcat(obj.Name, ' Smooth'))
            title(obj.Name)
            dynamicDateTicks
        end
        
        
        function plotRunningMeds(obj)
            hold on
            grid on
            plot(obj.DataObject.Time, obj.runningFSWmedian, 'y*')
            plot(obj.DataObject.Time, obj.runningTSWmedian, 'c*')
            legend(obj.Name, 'running FSW median', 'running TSW median');
            dynamicDateTicks
            hold off
            grid off
           
        end
        function plotInitFSWTransitions(obj)
             hold on
             grid on
             dynamicDateTicks;
             transitions = obj.getVar('name', 'TSW', 'level', ...
                 'raw', 'data', 'FSWtransitions');
            plot(transitions.StartTime, transitions.StartData, 'yd', ...
                'MarkerEdgeColor','k', 'MarkerFaceColor', 'yellow', 'MarkerSize',5);
            plot(transitions.EndTime, transitions.EndData, 'md', ...
                'MarkerEdgeColor','k', 'MarkerFaceColor','magenta', 'MarkerSize',5);
        end
        
        function plotInitTSWTransitions(obj)
             hold on
             grid on
             dynamicDateTicks;
             transitions = obj.getVar('name', 'TSW', 'level', ...
                 'raw', 'data', 'TSWtransitions');
             plot(transitions.StartTime, transitions.StartData, 'yo', ...
                'MarkerEdgeColor','k', 'MarkerFaceColor', 'yellow', 'MarkerSize',5);
            plot(transitions.EndTime, transitions.EndData, 'mo', ...
                'MarkerEdgeColor','k', 'MarkerFaceColor','magenta', 'MarkerSize',5);
        end
        
        function plotTransitionsVertical(obj)
             hold on
             grid on
             dynamicDateTicks;
            for i=1:length(obj.TransStartTime)
                line([obj.TransStartTime(i) obj.TransStartTime(i)], [0 2], 'color', 'green')
            end
            for i=1:length(obj.TransEndTime)
                line([obj.TransEndTime(i) obj.TransEndTime(i)], [0 2], 'color', 'red')
            end
        end
        
        function plotFSWSuspectTransitions(obj)
            plot(obj.TransStartTime(obj.TransStartFlag == 3), ...
                obj.TransStartData(obj.TransStartFlag == 3), ...
                'g*', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'green', ...
                'MarkerSize', 5);
        end
        
                
        function plotTSWSuspectTransitions(obj)
            plot(obj.TSWStartTime(obj.TSWStartFlag == 3), ...
                obj.TSWStartData(obj.TSWStartFlag == 3), ...
                'g*', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'green', ...
                'MarkerSize', 5);
        end
        
        function plotFSWGoodTransitions(obj)
            transitions = obj.getVar('name', 'TSW', 'level', ...
                 'preprocessed', 'data', 'FSWtransitions');
            plot(transitions.StartTime, transitions.StartData, ...
                'gd','MarkerEdgeColor','k','MarkerFaceColor', 'green', 'MarkerSize',5);
            plot(transitions.EndTime, transitions.EndData, ...
                 'cd', 'MarkerEdgeColor','k','MarkerFaceColor', 'cyan' ,'MarkerSize',5);
        end
        
        
        function plotTSWGoodTransitions(obj)
            transitions = obj.getVar('name', 'TSW', 'level', ...
                 'preprocessed', 'data', 'TSWtransitions');
            plot(transitions.StartTime, transitions.StartData, ...
                'go','MarkerEdgeColor','k','MarkerFaceColor','green','MarkerSize',5);
            plot(transitions.EndTime, transitions.EndData, ...
                'co','MarkerEdgeColor','k','MarkerFaceColor','cyan','MarkerSize',5);
        end
        
        
        function plotSyncData(obj)
            plot(obj.SyncDataObject.Time, obj.SyncDataObject.Data(:,20), 'k');
            legend(obj.Name);
        end
        

        
        
    end   % end methods
end   % end classDef