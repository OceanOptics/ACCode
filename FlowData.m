%FlowData - Data object to hold the all processing information for flow data
%FlowData inherits from AncillaryData
%
% Syntax:  flowdata = FlowData(nameIn, dataValuesIn, timestampsIn, unitsIn)
% Inputs:
%    nameIn - Description
%    dataValuesIn - Description
%    timestampsIn - Description
%    unitsIn - 
%
%
% Example: 
%    fd = FlowData('Flow1', Flow, fullDate, obj.FlowUnits);
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
classdef FlowData < AncillaryData
    properties
        Name
        Type = 'Flow';
        Units   % units of measurement for flow
        
        DataObject         % a timeseries object
        
        % preprocessing variables
        SmoothData         %smoothed data for finding transitions
        PPTimespan = 4800; %in seconds 4800=80 min
        SamplingFreq = 60; %per minute 60== 1/sec
        FSWFreq = 50;      %number of minutes 50 == every 50 minutes
        FSWDuration = 10;  %duration of filtered period in minutes
        TransStartData;    %start dataooint for FSW periods
        TransStartTime;    %start timestamp for FSW periods
        TransEndData;      %end datapoint for FSW periods
        TransEndTime;      %end timestamp for FSW periods
        TransStartFlag;    %flags for FSW start periods
        FSWPeriodLengths;  %lengths of each identified FSW period
        TransStartTimeGood;%checked as good FSW start periods
        TransEndTimeGood;  %checked as good FSW end periods
        levelsMap;
    end
    properties (Access = private)
        runningFSWmedian   % the running median used for finding transitions
        runningTSWmedian   % the running median used for finding transitions
    end    
    methods
        
        % constructor
        function obj = FlowData(nameIn, dataValuesIn, timestampsIn, unitsIn)
            
            %%% Pre-initialization %%%
            % Any code not using output argument (obj)
            
            %%% Object Initialization %%%
            % Call superclass constructor before accessing object
            % This statment cannot be conditionalized
                
            obj = obj@AncillaryData(nameIn, dataValuesIn, timestampsIn, unitsIn);
            
            %%% Post-initialization %%%
            % Any code, including access to the object
            % setInfo(obj, nameIn);
        end
        
        %%# checkTransitions
        function checkTransitions(obj)
             L = log4m.getLogger();
            %% check transitions are within expected time of each other
           
            % for each time in the list of start times
            uppermargin = datenum(0,0,0,0, obj.FSWFreq + .1*obj.FSWFreq, 0);
            lowermargin = datenum(0,0,0,0, obj.FSWFreq - .1*obj.FSWFreq, 0);

            if length(obj.TransStartTime) > 1
            
                   %check first one
                    if obj.TransStartTime(2) - obj.TransStartTime(1) < lowermargin || ...
                            obj.TransStartTime(end) - obj.TransStartTime(end-1) > uppermargin
                        obj.TransStartFlag(1) = 3; % suspect
                    else
                        obj.TransStartFlag(1) = 1; % good
                    end
                    sprintf('first: %u', obj.TransStartFlag(1))
                    % if there are more than 2
                    if length(obj.TransStartTime) >= 2
                        % for everyone after first one except last one
                        for i=2:(length(obj.TransStartTime) - 1)
                            sprintf('i: %u, start time: %s', i, datestr(obj.TransStartTime(i)))
                            % find the difference between this one and the previous
                            % one
                            timediffprev = obj.TransStartTime(i) - obj.TransStartTime(i-1);
                            timediffnext = obj.TransStartTime(i+1) - obj.TransStartTime(i);
                            sprintf('timediffprev: %s', datestr(timediffprev))
                            sprintf('timediffnext: %s', datestr(timediffnext))

                            if (timediffprev < uppermargin) && (timediffprev > lowermargin ) || ...
                                (timediffnext < uppermargin && timediffnext > lowermargin )
                                obj.TransStartFlag(i) = 1; %good
                            else
                                obj.TransStartFlag(i) = 3; %suspect
                            end
                            sprintf('flag: %u',  obj.TransStartFlag(i))
                            sprintf('----------------------------------')
                        end  % end for loop
                    end  % end if more than 2
                   %check last one
                   if obj.TransStartTime(end) - obj.TransStartTime(end-1) < lowermargin || ...
                            obj.TransStartTime(end) - obj.TransStartTime(end-1) > uppermargin
                        obj.TransStartFlag(end) = 3; % suspect
                    else
                        obj.TransStartFlag(end) = 1; % good
                   end
                   sprintf('last: %u', obj.TransStartFlag(end))

                   %% check transitions are length we expect
                   if obj.TransStartTime(1) < obj.TransEndTime(1)
                       disp('data starts in TSW, not FSW')

                        if length(obj.TransStartTime) == length(obj.TransEndTime)
                       %     disp('% we have complete set of transitions(?)')
                            obj.FSWPeriodLengths = (obj.TransEndTime - obj.TransStartTime)
                        elseif length(obj.TransStartTime) > length(obj.TransEndTime)
                        %    disp('we have an extra start time, i.e. data ends during')
                            % FSW

                            obj.FSWPeriodLengths = obj.TransEndTime - obj.TransStartTime(1:end-1);

                        else
                         %   disp('we shouldn''t have an extra end period when it starts in TSW')
                        end   % end length check
                   else
                       %obj.TransStartTime(1) > obj.TransEndTime(1)  % else start < end
                       %disp('% data starts in FSW')
                       if length(obj.TransStartTime) == length(obj.TransEndTime)
                        %   disp('it starts in fSW and ends in FSW')
                       elseif length(obj.TransStartTime) < length(obj.TransEndTime);
                         %  disp('% we end in TSW')
                       else
                          % disp('we shouldn''t have an extra start period')
                       end
                   end
                   % set margin
                     FSWuppermargin = datenum(0,0,0,0, obj.FSWDuration + .2*obj.FSWDuration, 0);
                     FSWlowermargin = datenum(0,0,0,0, obj.FSWDuration - .2*obj.FSWDuration, 0);
                     % check each duration is within margin
                     sprintf('lower limit %s', datestr(FSWuppermargin))
                     sprintf('upper limit %s', datestr(FSWlowermargin))
                   for i=1:length(obj.FSWPeriodLengths)
                       % if length of FSW is smaller than supposed duration -
                       % margin OR greater than supposed duration + margin
                       sprintf('period length: %s', datestr(obj.FSWPeriodLengths(i)))
                       if (obj.FSWPeriodLengths(i) < (FSWlowermargin) || ...
                               (obj.FSWPeriodLengths(i) > (FSWuppermargin)))
                           obj.TransStartFlag(i) = 3;
                       end
                   end
            elseif length(obj.TransStartTime) == 1 % length(obj.TransStartTime) > 1
                disp('only have one transitions times');
                obj.TransStartFlag(1) = 1;
            else
                L.error('FlowData', 'NO flow transitions times');
            end;
           
           
        end   % end function
        
        %%# setGoodTransitions
        function setGoodTransitions(obj)
            %assign transitions with a good flag
            obj.TransStartTimeGood = obj.TransStartTime(obj.TransStartFlag == 1);
        end
            
        %%# setTransitions
        function setTransitions( obj, TransStartDataIn, TransEndDataIn, TransStartTimeIn, TransEndTimeIn )
            %assign TSW transitions from external source
            obj.TransStartData = TransStartDataIn;
            obj.TransEndData = TransEndDataIn;
            obj.TransStartTime = TransStartTimeIn;
            obj.TransEndTime = TransEndTimeIn;
        
        end

        %%# setRunningMeds
        function setRunningMeds(obj)
            [obj.runningFSWmedian, obj.runningTSWmedian] = ...
                findFSWTSWRunMeds(obj.SmoothData, obj.DataObject.Time, obj.PPTimespan, obj.SamplingFreq);
        end
        
        %%# plotTransitions
        function plotTransitions(obj)
            
            hold on;
            grid on;
            plot(obj.DataObject.Time, obj.DataObject.Data, 'c')
            plot(obj.TransStartTime, obj.TransStartData, 'ks')
            plot(obj.TransEndTime, obj.TransEndData, 'kd')
            
            for i=1:length(obj.TransStartTime)
                line([obj.TransStartTime(i) obj.TransStartTime(i)], [0 35], 'color', 'green')
            end
            for i=1:length(obj.TransEndTime)
                line([obj.TransEndTime(i) obj.TransEndTime(i)], [0 35], 'color', 'red')
            end
            set(gca, 'YLim', [0 45])
            legend('flow and valve state', 'valve open', 'valve close')
            dynamicDateTicks
        end
        

        
        function plotSuspectTransitions(obj)
            plot(obj.TransStartTime(obj.TransStartFlag == 3), obj.TransStartData(obj.TransStartFlag == 3), 'k*')
        end
        
        function plotRunningMeds(obj)
            hold on
            grid on
            plot(obj.DataObject.Time, obj.SmoothData)
            plot(obj.DataObject.Time, obj.runningFSWmedian, 'y*')
            plot(obj.DataObject.Time, obj.runningTSWmedian, 'c*')
            legend(obj.Name, 'running FSW median', 'running TSW median');
            dynamicDateTicks
            hold off
            grid off
           
        end    
        
        function qaqc(obj)
            disp('Stub for qaqc():');
            obj.Name
        end   % end QAQC()
        
        function bin(obj)
            disp('Stub bin():');
            obj.Name
        end   % end Bin()
        
         %% find transitions from running medians          
        function findTransitions( obj ) 


            % create a logical array of when the TSWmedian has no data/is Nan
            % should create an array of 1's for when we have 'filtered' data
            nanArray = isnan(obj.runningTSWmedian);

            % get rid of minor fluctuations across the transition point
            % smooth across 15 seconds (4 measurements/second*15s = 60 measurements)
            % should create an array of 1s and 0s and values in between
            rm = smooth(nanArray,101);

            % get rid of rounding error
            rm_round = (floor(rm*100))/100;
            rm = rm_round;

            % create index of values where it is just 1/filtered, not transition
            idxFiltered = (rm == 1);

            diffArray =  idxFiltered(2:end) - idxFiltered(1:end-1);
            diffArray(end+1) = NaN;
            startsIndex = find(diffArray > 0)
            endsIndex = find(diffArray < 0)

            meanVecLength = length(obj.runningTSWmedian);
            startsLength = length(startsIndex);
            endsLength = length(endsIndex);

            startData = zeros( size( startsIndex));
            startData(:,:) = NaN;
            startTimes = zeros( size( startsIndex));
            startTimes(:,:) = NaN;

            endData = zeros( size( endsIndex));
            endData(:,:) = NaN;
            endTimes = zeros( size( endsIndex));
            endTimes(:,:) = NaN;

            % go through the starts and adjust for the lag created by 'smooth'
            for iStarts = 1:startsLength;
                thisStart = startsIndex(iStarts);

                if thisStart - 51 < 1
                    runningMedianIn(thisStart - 51)
                    startData(iStarts) = obj.runningTSWmedian(1);
                    startData(iStarts)
                    startTimes(iStarts) = obj.DataObject.Time(1);
                    datestr(startTimes(iStarts))
                else   % the start we'd want is past the end of the median vector
                    startData(iStarts) = obj.runningTSWmedian(thisStart -51);
                    startTimes(iStarts) = obj.DataObject.Time(thisStart-51);
                    datestr(startTimes(iStarts))
                end
            end

            for iEnds = 1:endsLength;
                % if the first end is before the beginning of the file (shouldn't
                % happen)
                thisEnd = endsIndex(iEnds);
                if (thisEnd + 51) < meanVecLength
                    endData(iEnds) = obj.runningTSWmedian(thisEnd+51);
                    endTimes(iEnds) = obj.DataObject.Time(thisEnd+51);
                    datestr(endTimes(iEnds))
                else
                    endData(iEnds) = obj.runningTSWmedian(end);
                    endTimes(iEnds) = obj.DataObject.Time(end);
                    datestr(endTimes(iEnds))
                end
            end
            obj.TransStartData = startData;
            obj.TransEndData = endData;
            obj.TransStartTime = startTimes;
            obj.TransEndTime = endTimes;
        end
        
        
        function plotData(obj)
            
            ts = obj.DataObject;
            plot(obj.DataObject.Time, obj.DataObject.Data)
            xlabel('Timestamps')
            ylabel('Flow')
            dynamicDateTicks();
            
        end   % end plotData
        
    end   % end methods
end   % end classDef