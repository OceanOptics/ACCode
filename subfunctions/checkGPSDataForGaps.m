function [binned_data_out, bin_flags_out] = checkGPSDataForGaps(unbinned_data, unbinned_timestamps, binned_data, binned_timestamps, bin_flags)
%checkGPSDataForGaps - Look for gaps in the GPS data and copy in data where
%appropriate
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [binned_data_out, binned_flags_out] = checkGPSDataForGaps(unbinned_data, unbinned_timestamps, binned_data, binned_timestamps, bin_flags)
%
% Inputs:
%    unbinned_data - raw data
%    unbinned_timestamps - raw timestamps
%    binned_data - data binned to fixed set of timestamps
%    binned_timestamps - timestamps of bins
%    bin_flags - any existing flags for the bins
%
% Outputs:
%    binned_data_out - bins with data changes made
%    bin_flags_out - the accompanying flags for data changed
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TemperatureData,  GPSData, SalinityData

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% Nov 2015; Last revision: 13-Nov-15

%------------- BEGIN CODE --------------
    
    L = log4m.getLogger();

    time = unbinned_timestamps;
    data = unbinned_data;
    bins = binned_timestamps;
    bin_data = binned_data;
    flags = bin_flags;

    params.GPS_GAP_THRESHOLD = 2;
    lowThreshold = params.GPS_GAP_THRESHOLD;
 

          
    % produce datevec (YMDHMS) of the diff between each time and
    % next
    timestampDiff = datevec(diff(time));   %obj.DataObject.Time

    % sum hours and minutes to create total minutes of difference
    totalMinutesDiff = timestampDiff(:,4)*.60 + timestampDiff(:,5);

    % create index to timestamps with difference greater
    % than 2 minutes to the next one
    minuteDiffIndex = totalMinutesDiff >= lowThreshold;

    % check if any timestamps have big difference
    if sum(minuteDiffIndex) >= 1
        
        L.debug('checkGPSDataForGaps', 'Gap in data, flagging data');
        
%                 % set the quality flag for this datapoint to 3/suspect
%                 allData.TemperatureData.DataObject.Quality(minuteDiffIndex) = 3;
    else
%                 obj.L.debug('TemperatureData.bin()', 'No large gaps in data');
        L.debug('checkGPSDataForGaps', 'no large gaps in data');
    end;
            
%             % 2.  Mark bins with suspect datapoints
% 
%             % create the bin Flags; Set to '1' by default ('good')
% %             obj.BinFlags = ones(size(obj.BinnedTimestamps));

    % find the timestamps which come before and after the gap in data
    tsBeforeGapIndex = find(minuteDiffIndex);
    tsBeforeGapTimestamps = time(minuteDiffIndex); 
    tsAfterGapIndex = find(minuteDiffIndex) + 1;
    tsAfterGapTimestamps = time(tsAfterGapIndex);

%     datestr(tsBeforeGapTimestamps)
%     datestr(tsAfterGapTimestamps)
      %%      
%     disp('****************************************************')
    for iTimestamp = 1:size(tsBeforeGapTimestamps)
%         disp('===================================================')
%         iTimestamp
%         disp('===================================================')

        thisTimestamp = tsBeforeGapTimestamps(iTimestamp);
        thisTSIndex = tsBeforeGapIndex(iTimestamp);
        nextTimestamp = tsAfterGapTimestamps(iTimestamp);
        nextTSIndex = tsAfterGapIndex(iTimestamp);

        % check this and next TS are within bins we're using
        if thisTimestamp >= bins(1) && thisTimestamp <= bins(end)

%             disp('thisTS')
%             datestr(thisTimestamp)
%             disp('data')
%             data(thisTSIndex)
%             disp('nextTS')
%             datestr(nextTimestamp)
%             data(nextTSIndex)
%             disp('------------')
% 
%             disp('time between')
            timestampDiff = datevec(nextTimestamp - thisTimestamp);

            % sum hours and minutes to create total minutes of difference
            totalMinutesDiff = timestampDiff(:,4)*60 + timestampDiff(:,5);
%             totalMinutesDiff

            % if total minutes in gap is 120 or less; split gap in two
            if totalMinutesDiff <= 120
%                 disp('gap < 120')
% 
%                 disp('half gap is ')
%                 totalMinutesDiff/2
%                 disp('midway is')
                midTimestamp = addtodate(thisTimestamp, floor(totalMinutesDiff/2), 'minute');
%                 datestr(midTimestamp)

                % for first half, copy in first data point

                indexBin = bins > thisTimestamp & bins < midTimestamp;
%                 binNumbers = find(indexBin)
%                 bin_data(binNumbers)
                bin_data(indexBin) = data(thisTSIndex);
                bin_flags(indexBin) = 4;  % interpolated
%                 bin_data(binNumbers)              

                % for second half, copy in last data point
                indexBin = bins > midTimestamp & bins < nextTimestamp;

                binNumbers = find(indexBin);
%                 bin_data(binNumbers)                   
                bin_data(indexBin) = data(nextTSIndex);
                bin_flags(indexBin) = 4;  % interpolated               
%                 bin_data(binNumbers)  

            else
%                 disp('gap > 120')
%                 disp('first end point is 1 hr after')
                firstHourTimestamp = addtodate(thisTimestamp, 60, 'minute');
%                 datestr(firstHourTimestamp)
%                 disp('second start point is 1 hr before')
                secondHourTimestamp = addtodate(nextTimestamp, -60, 'minute');
%                 datestr(secondHourTimestamp)


                % for first hour, copy in first data point
                indexBin = bins > thisTimestamp & bins < firstHourTimestamp;
                if sum(indexBin) > 1
%                     disp('found bin to change')
                    binNumbers = find(indexBin);
%                     bin_data(binNumbers)
                    bin_data(indexBin) = data(thisTSIndex);
                    bin_flags(indexBin) = 4;  % interpolated
%                     bin_data(binNumbers)
                else
%                     disp('no bins to find between')
%                     disp(datestr(thisTimestamp))
%                     disp('and')
%                     disp(datestr(firstHourTimestamp))
% 
%                     disp('first bin:')
%                     disp(datestr(bins(1)))
%                     disp('last bin:')
%                     disp(datestr(bins(end)))

                end;

                % for last hour, copy in last data point
                indexBin = bins > secondHourTimestamp & bins < nextTimestamp;
                if sum(indexBin) > 1
%                     disp('found bin to change')
%                     binNumbers = find(indexBin)
%                     bin_data(binNumbers)
                    bin_data(indexBin) = data(nextTSIndex);
                    bin_flags(indexBin) = 4;  % interpolated
%                     bin_data(binNumbers)
                else
%                     disp('no bins to find between')
%                     disp(datestr(secondHourTimestamp))
%                     disp('and')
%                     disp(datestr(nextTimestamp))
% 
%                     disp('first bin:')
%                     disp(datestr(bins(1)))
%                     disp('last bin:')
%                     disp(datestr(bins(end)))

                end;
            end

        elseif thisTimestamp < bins(1) && nextTimestamp < bins(1)
%             disp('both timestamps in this gap come before bins start')
%             datestr(bins(1))
%             disp('and end')
%             datestr(bins(end))
        elseif thisTimestamp > bins(end)
%                 disp('first timestamp in this gap comes after bins end')
%                 datestr(bins(end))
        else
%             disp('something wierd going on with timestamps here')
        end;

    end
    
    binned_data_out = bin_data;
    bin_flags_out = bin_flags;

end

%------------- END --------------


%%
            
% %                 indexSuspectBin = obj.var.(level).BinnedTimestamps < thisTimestamp &...
% %                     obj.var.(level).BinnedTimestamps > prevmin;
%                 indexSuspectBin = bins < thisTimestamp & bins > nextTimestamp;
% 
%                 if sum(indexSuspectBin) >=1
% %                     obj.L.debug('TemperatureData.bin()', 'Found suspect bins')
% disp('suspect bins found')
% sum(indexSuspectBin)
% % 
% %                     obj.var.(level).LabTempBinFlags(indexSuspectBin) = 3;
% %                     obj.var.(level).BoatTempBinFlags(indexSuspectBin) = 3;
%                 else
%                     obj.L.debug('TemperatureData.bin()', 'No suspect bins found')
%                     % do nothing
%                 end;
%             end;