%%ACData Class

classdef ACData < AncillaryData
    properties
        Name
        Type = 'ac';
        Units
        
        % a timeseries object
        DataObject
        SyncDataObject
        
        % preprocessing variables
        SmoothData
        PPTimespan = 19200;
        SamplingFreq = 240;
        
        FSWFreq = 60;      %number of minutes 50 == every 50 minutes
        TSWDuration = 50;
        FSWDuration = 10;
        
        
        TransStartData;
        TransStartTime;
        TransEndData;
        TransEndTime;
        TransStartFlag;
        TransEndFlag;
        FSWPeriodLengths;
        
        GoodFSWStartTimes;
        GoodFSWStartData;
        GoodFSWEndTimes;
        GoodFSWEndData;
    
        
        TSWStartData;
        TSWStartTime;
        TSWEndData;
        TSWEndTime;
        TSWStartFlag;
        TSWEndFlag;
        TSWPeriodLengths;
                
        GoodTSWStartTimes;
        GoodTSWStartData;
        GoodTSWEndTimes;
        GoodTSWEndData;
        
        TSWIndex;
        FSWIndex;
        levelsMap;
        
        L
        runningFSWmedian
        runningTSWmedian
        
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
                
            obj = obj@AncillaryData(nameIn, dataValuesIn, timestampsIn, unitsIn);
            
            %%% Post-initialization %%%
            % Any code, including access to the object
            % adding variables to timeseries object.
            obj.L = log4m.getLogger();
            obj.L.debug('ACData.ACData()','Created object');
            %setInfo(obj, nameIn);
        end
        %% ----------------------------------------------------------------
        % setGoodTransitions
        % -----------------------------------------------------------------
        function setGoodFSWTransitions(obj, goodStarts, goodEnds)
            
            goodStarts = goodStarts(~isnan(goodStarts));
            goodEnds = goodEnds(~isnan(goodEnds));
               
            % now only use good start times
%             size(obj.TransStartTime)
%             size(goodStarts)
%             size(obj.TransEndTime)
%             size(goodEnds)
            
%             % create an index to the good times/data in the original list
%             [~,ind1] = min(abs(bsxfun(@minus,datenum(obj.TransStartTime),datenum(goodStarts)')));
%             % subset the data to the good times
%             goodStartData = obj.TransStartData(ind1,:);

            [~,ind1] = min(abs(bsxfun(@minus,datenum(obj.DataObject.Time),datenum(goodStarts)')));
            % subset the data to the good times
            goodStartData = obj.DataObject.Data(ind1,20); 

            % have goodStarts, now need good 

            
            % do same for end times
%             [~,ind1] = min(abs(bsxfun(@minus,datenum(obj.TransEndTime),datenum(goodEnds)')));
%             goodEndData = obj.TransEndData(ind1,:);

            [~,ind1] = min(abs(bsxfun(@minus,datenum(obj.DataObject.Time),datenum(goodEnds)')));
            goodEndData = obj.DataObject.Data(ind1,20);
%             size(ind1)
%             size(goodEndData)
            
            
            
            
            
            % set the object properties manually
            obj.GoodFSWStartData = goodStartData;
            obj.GoodFSWStartTimes = goodStarts;
            obj.GoodFSWEndData = goodEndData;
            obj.GoodFSWEndTimes = goodEnds;
             
%             for i = 1:length(obj.GoodFSWStartTimes)
%                 sprintf('%s %s', datestr(obj.GoodFSWStartTimes(i)), datestr(obj.GoodFSWEndTimes(i)))
%             end;
            
           
        end
        
        function setGoodTSWTransitions(obj, goodStarts, goodEnds)
            
            goodStarts = goodStarts(~isnan(goodStarts));
            goodEnds = goodEnds(~isnan(goodEnds));
               
            % now only use good start times
%             size(obj.TSWStartTime)
%             size(goodStarts)
%             size(obj.TSWEndTime)
%             size(goodEnds)
            
            % find matching data to timestamps:
            
            % create an index to the good times/data in the original list
%             [~,ind1] = min(abs(bsxfun(@minus,datenum(obj.TSWStartTime),datenum(goodStarts)')));
%             % subset the data to the good times
%             goodStartData = obj.TSWStartData(ind1,:);
            [~,ind1] = min(abs(bsxfun(@minus,datenum(obj.DataObject.Time),datenum(goodStarts)')));
            % subset the data to the good times
            goodStartData = obj.DataObject.Data(ind1,20); 


            % do same for end times
%             [~,ind1] = min(abs(bsxfun(@minus,datenum(obj.TSWEndTime),datenum(goodEnds)')));
%             goodEndData = obj.TSWEndData(ind1,:);
            [~,ind1] = min(abs(bsxfun(@minus,datenum(obj.DataObject.Time),datenum(goodEnds)')));
            goodEndData = obj.DataObject.Data(ind1,20);
            
            % set the object properties manually
            obj.GoodTSWStartData = goodStartData;
            obj.GoodTSWStartTimes = goodStarts;
            obj.GoodTSWEndData = goodEndData;
            obj.GoodTSWEndTimes = goodEnds;
             
%              for i = 1:length(obj.GoodTSWStartTimes)
%                 sprintf('%s %s', datestr(obj.GoodTSWStartTimes(i)), datestr(obj.GoodTSWEndTimes(i)))
%             end;
            
           
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
            obj.L.info('ACData.setSmoothData', 'start method');
        end

        %% ----------------------------------------------------------------
        % setRunningMedians
        % calculate the running medians - used to identify TSW and FSW
        % -----------------------------------------------------------------
        function setRunningMeds(obj)
            obj.L.info('ACData.setRunningMeds', 'start method');
            [obj.runningFSWmedian, obj.runningTSWmedian] = ...
                findFSWTSWRunMeds(obj.SmoothData(:,20), obj.DataObject.Time, obj.PPTimespan, obj.SamplingFreq);
            obj.L.info('ACData.setRunningMeds', 'end method');
        end
        
        
        %% ----------------------------------------------------------------
        % findTransitions
        % this is for filtered FSW
        % -----------------------------------------------------------------
        function findTransitions(obj)
            
            obj.L.info('ACData.findTransitions', 'start method');
            
            nanArray = isnan(obj.runningTSWmedian);
            obj.L.debug('ACData.findTransitions', ...
                sprintf('sum(nanArray): %u', sum(nanArray)));
            obj.L.debug('ACData.findTransitions', ...
                sprintf('size of nanArray: %u x %u', size(nanArray)));

            % get rid of minor fluctuations across the transition point
            % smooth across 15 seconds (4 measurements/second*15s = 60 measurements)
            % should create an array of 1s and 0s and values in between
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
            
            % same to here
            obj.L.debug('ACData.findFSWTransitions',...
                sprintf('length of running median: %u', length(obj.runningFSWmedian)));
            obj.L.debug('ACData.findFSWTransitions',...
                sprintf('length of startsIndex: %u', length(startsIndex)));
            for i = 1:length(startsIndex)
                obj.L.debug('ACData.findFSWTransitions',...
                    sprintf('startsIndex(i)+76: %u', startsIndex(i) + 76));
                if startsIndex(i) + 76 > length(obj.runningFSWmedian)
                    obj.L.error('ACData.findFSWTransitions',...
                        sprintf('About to run over end of runningFSWmedian; i: %u', i));
                    % make this the end
                    obj.TransStartData(i) = obj.runningFSWmedian(startsIndex(end));
                    obj.TransStartTime(i) = obj.DataObject.Time(startsIndex(end));
                   
                else
                    obj.L.debug('ACData.findFSWTransitions','Starts index + 76 within end');
    %             obj.FSWStartData = obj.runningFSWmedian(startsIndex + 76);
    %             obj.FSWStartTime = obj.DataObject.Time(startsIndex + 76);
    %             obj.FSWEndData = obj.runningFSWmedian(endsIndex-76);
    %             obj.FSWEndTime = obj.DataObject.Time(endsIndex-76);
                    obj.TransStartData(i) = obj.runningFSWmedian(startsIndex(i) + 76);
                    obj.TransStartTime(i) = obj.DataObject.Time(startsIndex(i) + 76);

                end
            end  % end startsIndex loop
            for i = 1:length(endsIndex)
                obj.L.debug('ACData.findFSWTransitions',...
                    sprintf('endsIndex(i)-76: %u', endsIndex(i) - 76));
                if endsIndex(i) - 76 < 1
                    obj.L.error('ACData.findFSWTransitions',...
                        sprintf('About to run over beginning of runningFSWmedian; i: %u', i));
                    % make this beginning 
                     obj.TransEndData(i) = obj.runningFSWmedian(endsIndex(1));
                     obj.TransEndTime(i) = obj.DataObject.Time(endsIndex(1));
                else
                     obj.L.debug('ACData.findFSWTransitions','ends index - 76 within beginning'); 
                     obj.TransEndData(i) = obj.runningFSWmedian(endsIndex(i)-76);
                     obj.TransEndTime(i) = obj.DataObject.Time(endsIndex(i)-76);
                end
            end  % end endsIndex loop
        
            % how did these get turned around????
            obj.TransStartTime = obj.TransStartTime';
            obj.TransEndTime = obj.TransEndTime';
            obj.TransStartData = obj.TransStartData';
            obj.TransEndData = obj.TransEndData';
            
            % just being set to 1
            obj.TransStartFlag = ones(size(obj.TransStartTime));
            obj.TransEndFlag = ones(size(obj.TransEndTime));
        end
        
        %% ----------------------------------------------------------------
        % findTSWTransitions
        % -----------------------------------------------------------------
        function findTSWTransitions(obj)
            
            obj.L.info('ACData.findTSWTransitions', 'start method');
            
            nanArray = isnan(obj.runningFSWmedian);

            % get rid of minor fluctuations across the transition point
            % smooth across 15 seconds (4 measurements/second*15s = 60 measurements)
            % should create an array of 1s and 0s and values in between
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
            
            obj.L.debug('ACData.findTSWTransitions',...
                sprintf('length of running median: %u', length(obj.runningTSWmedian)));
            obj.L.debug('ACData.findTSWTransitions',...
                sprintf('length of startsIndex: %u', length(startsIndex)));
            for i = 1:length(startsIndex)
                obj.L.debug('ACData.findTSWTransitions',...
                    sprintf('startsIndex(i)+76: %u', startsIndex(i) + 76));
                if startsIndex(i) + 76 > length(obj.runningTSWmedian)
                    obj.L.error('ACData.findTSWTransitions',...
                        sprintf('About to run over end of runningTSWmedian; i: %u', i));
                    % make this the end
                    obj.TSWStartData(i) = obj.runningTSWmedian(startsIndex(end));
                    obj.TSWStartTime(i) = obj.DataObject.Time(startsIndex(end));
                   
                else
                    obj.L.debug('ACData.findTSWTransitions','Starts index + 76 within end');
    %             obj.TSWStartData = obj.runningTSWmedian(startsIndex + 76);
    %             obj.TSWStartTime = obj.DataObject.Time(startsIndex + 76);
    %             obj.TSWEndData = obj.runningTSWmedian(endsIndex-76);
    %             obj.TSWEndTime = obj.DataObject.Time(endsIndex-76);
                    obj.TSWStartData(i) = obj.runningTSWmedian(startsIndex(i) + 76);
                    obj.TSWStartTime(i) = obj.DataObject.Time(startsIndex(i) + 76);

                end
            end  % end startsIndex loop
            for i = 1:length(endsIndex)
                obj.L.debug('ACData.findTSWTransitions',...
                    sprintf('endsIndex(i)-76: %u', endsIndex(i) - 76));
                if endsIndex(i) - 76 < 1
                    obj.L.error('ACData.findTSWTransitions',...
                        sprintf('About to run over beginning of runningTSWmedian; i: %u', i));
                    % make this beginning 
                     obj.TSWEndData(i) = obj.runningTSWmedian(endsIndex(1));
                     obj.TSWEndTime(i) = obj.DataObject.Time(endsIndex(1));
                else
                    obj.L.debug('ACData.findTSWTransitions','ends index - 76 within beginning'); 
                     obj.TSWEndData(i) = obj.runningTSWmedian(endsIndex(i)-76);
                    obj.TSWEndTime(i) = obj.DataObject.Time(endsIndex(i)-76);
                end
            end  % end endsIndex loop
        
            % how did these get turned around????
            obj.TSWStartTime = obj.TSWStartTime';
            obj.TSWEndTime = obj.TSWEndTime';
            obj.TSWStartData = obj.TSWStartData';
            obj.TSWEndData = obj.TSWEndData';
            
            % just being set to one
            obj.TSWStartFlag = ones(size(obj.TSWStartTime));
            obj.TSWEndFlag = ones(size(obj.TSWEndTime));
        end    % end findTSWTransitions

        %% ---------------------------------------------------------------
        % ----------------------------------------------------------------
        % Sync to should be run AFTER checkTransitions()
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        function syncTo(obj, AncillaryDataIn)
            % check transitions here
            % WE ARE MATCHING FILTERED STARTS HERE, NOT TOTAL
            goodTrans = obj.GoodFSWStartTimes;
            obj.L.debug('ACData.syncTo', ...
                sprintf('first good transition in FSW: %s', datestr(goodTrans(1))));
            
            if length(goodTrans) >= 1 
                obj.L.debug('ACData.syncTo', 'at least one goodTrans');
                dataInGood = AncillaryDataIn.TransStartTime(AncillaryDataIn.TransStartFlag == 1);
                dataInGood
                [~,ind1] = min(abs(datenum(dataInGood) - datenum(goodTrans(1))));
                ind1
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
                
%                 disp('1st good starts before')
%                 datestr(obj.GoodTSWStartTimes(1))
                
                
                obj.GoodTSWStartTimes = obj.GoodTSWStartTimes + offset;
                obj.GoodTSWEndTimes = obj.GoodTSWEndTimes + offset;
                obj.GoodFSWStartTimes = obj.GoodFSWStartTimes + offset;
                obj.GoodFSWEndTimes = obj.GoodFSWEndTimes + offset;
                
                
                
%                 disp('1st good starts after')
%                 datestr(obj.GoodTSWStartTimes(1))
                
            else
                obj.L.error('ACData.syncTo', 'No good transitions to sync to!');
            end
        end
               
        %% ----------------------------------------------------------------
        % -----------------------------------------------------------------
        % create two indexes TSWIndex, FSWIndex to mark TSW and FSW and
        % then flag transition periods as 'suspect' (3)
        function separateTSWFSW(obj)
            % create two NaN arrays the size of data
            
            % CHANGE THESE DATAOBJECTS BACKINTO DATASYNCOBJECTS
            % OR SUB DATASYNC FOR DATA
            
            %obj.TSWIndex = zeros( size( obj.SyncDataObject.Time ));
            % want 0s for false, not nans
            obj.FSWIndex = zeros( length( obj.SyncDataObject.Time ), 1);
            obj.TSWIndex = zeros( length( obj.SyncDataObject.Time ), 1);

            % get index of start times in data
            [~,indFSWStarts] = min(abs(bsxfun(@minus,datenum(obj.SyncDataObject.Time),datenum(obj.GoodFSWStartTimes)')));
            [~,indFSWEnds] = min(abs(bsxfun(@minus,datenum(obj.SyncDataObject.Time),datenum(obj.GoodFSWEndTimes)')));
             
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
            [~,indTSWStarts] = min(abs(bsxfun(@minus,datenum(obj.SyncDataObject.Time),datenum(obj.GoodTSWStartTimes)')));
            [~,indTSWEnds] = min(abs(bsxfun(@minus,datenum(obj.SyncDataObject.Time),datenum(obj.GoodTSWEndTimes)')));
             
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
        

        

        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        % offsetFSWStartTimes() 
        % This function will shift forward the FSW start times by a given
        % time (in seconds)
        function offsetFSWStartTimes(obj, offsetTimeIn)
        
            % NEED CHECK HERE FOR SYNCDATAOBJECT
            
            % convert time in seconds to a date vec
            offsetDatenum = datenum(0, 0, 0, 0, 0, offsetTimeIn);
            % for each start time:
            for iStartTime = 1:numel(obj.GoodFSWStartTimes)
            %   add offsetTime 
                obj.L.debug('ACData.offsetFSWStartTimes', ...
                    sprintf('old time: %s', datestr(obj.GoodFSWStartTimes(iStartTime))));
                
                tempStartTime = obj.GoodFSWStartTimes(iStartTime) + offsetDatenum;
                
                obj.L.debug('ACData.offsetFSWStartTimes', ...
                    sprintf('temp start time: %s', datestr(tempStartTime)));  
                
                %   find index of new start time
%               [~,ind1] = min(abs(datenum(dataInGood) - datenum(goodTrans(1) )));
                [~,ind1] = min(abs(datenum(obj.SyncDataObject.Time) - datenum( tempStartTime )));

                %   reasssign startData with correct time and data
                obj.GoodFSWStartTimes(iStartTime) = obj.SyncDataObject.Time(ind1,:); 
                
                % just use col 20 here -- that's what orig running med used
                obj.GoodFSWStartData(iStartTime) = obj.SyncDataObject.Data(ind1,20);

            end   % end for loop
            
        end  % end function offsetFSWStartTimes
        

       
        %% check calculated transitions
        % - rewrite this to just copy over transitions
        function checkTransitions(obj)
            
                    % loop through starts and ends
            startsLength = length(obj.TransStartTime);
            endsLength = length(obj.TransEndTime);
            
            % which is longer?
            if startsLength >= endsLength
                maxLength = startsLength + 1;
            else
                maxLength = endsLength + 1;
            end;
            
            % set up flags
            obj.TransStartFlag = zeros( startsLength, 1 );
            obj.TransStartFlag(:,:) = NaN;
            
            % set up two new start/end arrays
            tempEnds = zeros(maxLength, 1);
            tempEnds(:,:) = NaN;  
            
            % -- COPIED FROM CHECK TRANSTIONS ABOVE --
            % check how many ends are before the first beginning 
            % i.e. start in middle of period.
            % if yes, make first start a NAN
            tempidx = obj.TransEndTime < obj.TransStartTime(1);
            numEndsWithNoStarts = sum(tempidx);
            obj.L.debug('ACData.checkFSWTransitions', ...      
               sprintf('numEndsWithNoStarts: %u', numEndsWithNoStarts));
           
           % set up temparray of starts, padded by any nans needed to make
           % up for early ends before the first start
           if numEndsWithNoStarts > 0
               % if we have ends with no starts, increase the array of
               % starts by this number of ends
               tempStarts = zeros(numEndsWithNoStarts+length(obj.TransStartTime), 1);
               tempStarts(:,:) = NaN;
               obj.L.debug('ACData.checkFSWTransitions', ... 
                   sprintf('size of tempStarts: %u', length(tempStarts)));
               obj.L.debug('ACData.checkFSWTransitions', ... 
                   sprintf('length of start times: %u', length(obj.TransStartTime)));
               startAt = numEndsWithNoStarts+1;
               obj.L.debug('ACData.checkFSWTransitions', ...
                   sprintf('numEnds+1: %u',  startAt ));
               
               obj.L.debug('ACData.checkFSWTransitions', ...
                   sprintf('size(tempStarts): %u',  length(tempStarts(startAt:end)) ));               
               % we have some end times before our first start time
               tempStarts(startAt:end) = obj.TransStartTime(:,:);
               
               % do same for flags
               tempFlags = zeros(numEndsWithNoStarts+length(obj.TransStartTime), 1);
               tempFlags(:,:) = NaN;
               obj.L.debug('ACData.checkFSWTransitions', ...
                   sprintf('size numEnds+1: %u', length(tempFlags(numEndsWithNoStarts+1:end))));
               tempFlags(numEndsWithNoStarts+1:end) = obj.TransStartFlag;
               obj.L.debug('ACData.checkFSWTransitions', ...
                   sprintf('number good flags: %u', sum(tempFlags==1)));
               %mark the nan flags as good
               tempFlags(1:numEndsWithNoStarts) = 1;
                obj.L.debug('ACData.checkFSWTransitions', ...
                   sprintf('number good flags after fixing nans: %u', sum(tempFlags==1)));
               
               goodStarts = zeros(numEndsWithNoStarts+length(obj.TransStartTime), 1);
               goodStarts(:,:) = NaN;
               goodEnds = zeros(numEndsWithNoStarts+length(obj.TransStartTime), 1);
               goodEnds(:,:) = NaN; 
           else
               % no ends before first start time
               tempStarts = zeros(maxLength, 1);
               tempStarts(:,:) = NaN;
               tempStarts = obj.TransStartTime;
               tempFlags = zeros(maxLength, 1);
               tempFlags(:,:) = NaN;
               tempFlags = obj.TransStartFlag; 
               goodStarts = zeros(maxLength, 1);
               goodStarts(:,:) = NaN;
               goodEnds = zeros(maxLength, 1);
               goodEnds(:,:) = NaN; 
               
           end;
               
           % for each good flagged start time, find the corresponding end
           % time and copy it into a new array (don't worry about it
           % changing size on each iteration)
           % look for end marker at end of period
           obj.L.debug('ACData.checkFSWTransitions', ...
               sprintf('# tempStarts %u', length(tempStarts)));
           obj.L.debug('ACData.checkFSWTransitions', ...
                sprintf('# good flags %u', sum(tempFlags==1)));
            
           % only carry good starts across
           tempStarts = tempStarts(tempFlags == 1);
           obj.L.debug('ACData.checkFSWTransitions', ...
               sprintf('size of new tempStarts with only good data: %u', length(tempStarts)));
           
           for iStart = 1:size(tempStarts)
                obj.L.debug('ACData.checkFSWTransitions', ...
                    sprintf('looping through tempStarts: %u', iStart));
                
               if ~isnan(tempStarts(iStart))
                  obj.L.debug('ACData.checkFSWTransitions', ...
                      ('current one not nan'));
                   % find index of ends that fall within search area
                   % search area is this start +/- duration
                   aboveIndex = (obj.TransEndTime > (tempStarts(iStart) - lowerDurLimit));
                   belowIndex = (obj.TransEndTime < (tempStarts(iStart) + upperDurLimit));
                   searchAreaIndex = aboveIndex & belowIndex;
                   numEndsFound = sum(searchAreaIndex);
                  obj.L.debug('ACData.checkFSWTransitions', ...
                      sprintf('found %u ends in search area', numEndsFound));
                  
                   if numEndsFound == 1
                       % we have one exact match
                      obj.L.debug('ACData.checkFSWTransitions', ...
                          'we have exactly one match');
                       goodStarts(iStart) = tempStarts(iStart);
                       goodEnds(iStart) = obj.TransEndTime(searchAreaIndex);
                       obj.L.debug('ACData.checkFSWTransitions', ...
                           sprintf('copied over: start: %s; end: %s', datestr(goodStarts(iStart)), datestr(goodEnds(iStart))));
                   elseif numEndsFound == 0
                       % we have more than one match?
                       %
                      obj.L.debug('ACData.checkFSWTransitions', ...
                          'we have no matching end');
                      % don't copy across
                   else %numendsfound > 1
                       % we have more than one match
                       obj.L.debug('ACData.checkFSWTransitions', ...
                           'we have more than one match');
                       
                   end
               else  %is nan
                   goodStarts(iStart) = tempStarts(iStart);
                   goodEnds(iStart) = tempEnds(iStart);
               end
           end   % for loop iStarts
                   
           % get rid of nans
           goodStarts = goodStarts(~isnan(goodStarts));
           goodEnds = goodEnds(~isnan(goodEnds));
               
            % now only use good start times
            size(obj.TransStartTime)
            size(goodStarts)
            size(obj.TransEndTime)
            size(goodEnds)
             %[~,ind1] = min(abs(bsxfun(@minus,datenum(obj.TransStartTime),datenum(goodEnds)')));
             [~,ind1] = min(abs(bsxfun(@minus,datenum(obj.TransStartTime),datenum(goodStarts)')));
             goodStartData = obj.TransStartData(ind1,:);
             [~,ind1] = min(abs(bsxfun(@minus,datenum(obj.TransEndTime),datenum(goodEnds)')));
             goodEndData = obj.TransEndData(ind1,:);
             
             obj.GoodFSWStartData = goodStartData;
             obj.GoodFSWStartTimes = goodStarts;
             obj.GoodFSWEndData = goodEndData;
             obj.GoodFSWEndTimes = goodEnds;
             
             for i = 1:length(obj.GoodFSWStartTimes)
                sprintf('%s %s', datestr(obj.GoodFSWStartTimes(i)), datestr(obj.GoodFSWEndTimes(i)))
            end;
           
        end   % end function
        
        %%
        function qaqc(obj)
            if obj.Debug == 1;
                disp('Stub for qaqc():');
                obj.Name
            end;
        end   % end QAQC()
        
        function bin(obj)
            if obj.Debug == 1;
                disp('Stub bin():');
                obj.Name
            end;
        end   % end Bin()
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
%             plot(obj.DataObject.Time, obj.SmoothData(:,20))
            plot(obj.DataObject.Time, obj.runningFSWmedian, 'y*')
            plot(obj.DataObject.Time, obj.runningTSWmedian, 'c*')
            legend(obj.Name, 'running FSW median', 'running TSW median');
            dynamicDateTicks
            hold off
            grid off
           
        end
        function plotInitFSWTransitions(obj)
%             plot( obj.DataObject.Time, obj.DataObject.Data(:,20) )
             hold on
             grid on
             dynamicDateTicks;
            plot(obj.TransStartTime, obj.TransStartData, 'yd', ...
                'MarkerEdgeColor','k', 'MarkerFaceColor', 'yellow', 'MarkerSize',5);
            plot(obj.TransEndTime, obj.TransEndData, 'md', ...
                'MarkerEdgeColor','k', 'MarkerFaceColor','magenta', 'MarkerSize',5);
        end
        
        function plotInitTSWTransitions(obj)
%             plot( obj.DataObject.Time, obj.DataObject.Data(:,20) )
             hold on
             grid on
             dynamicDateTicks;
            plot(obj.TSWStartTime, obj.TSWStartData, 'yo', ...
                'MarkerEdgeColor','k', 'MarkerFaceColor', 'yellow', 'MarkerSize',5);
            plot(obj.TSWEndTime, obj.TSWEndData, 'mo', ...
                'MarkerEdgeColor','k', 'MarkerFaceColor','magenta', 'MarkerSize',5);
%             legend(obj.Name, 'starts', 'ends')
        end
        
        function plotTransitionsVertical(obj)
%             plot(obj.DataObject.Time, obj.DataObject.Data(:,20))
             hold on
             grid on
             dynamicDateTicks;
            for i=1:length(obj.TransStartTime)
                line([obj.TransStartTime(i) obj.TransStartTime(i)], [0 2], 'color', 'green')
            end
            for i=1:length(obj.TransEndTime)
                line([obj.TransEndTime(i) obj.TransEndTime(i)], [0 2], 'color', 'red')
            end
%             xlabel('timestamp')
%             ylabel(obj.Name)
%             legend(obj.Name, 'FSW start' , 'FSW end')
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
            plot(obj.GoodFSWStartTimes, obj.GoodFSWStartData, ...
                'gd','MarkerEdgeColor','k','MarkerFaceColor', 'green', 'MarkerSize',5);
            plot(obj.GoodFSWEndTimes, obj.GoodFSWEndData, ...
                 'cd', 'MarkerEdgeColor','k','MarkerFaceColor', 'cyan' ,'MarkerSize',5);
%              legend('good start', 'good end')
        end
        
        
        function plotTSWGoodTransitions(obj)
            plot(obj.GoodTSWStartTimes, obj.GoodTSWStartData, ...
                'go','MarkerEdgeColor','k','MarkerFaceColor','green','MarkerSize',5);
            plot(obj.GoodTSWEndTimes, obj.GoodTSWEndData, ...
                'co','MarkerEdgeColor','k','MarkerFaceColor','cyan','MarkerSize',5);
        end
        
        
        function plotSyncData(obj)
            plot(obj.SyncDataObject.Time, obj.SyncDataObject.Data(:,20), 'k');
            legend(obj.Name);
        end
        

        
        
    end   % end methods
end   % end classDef