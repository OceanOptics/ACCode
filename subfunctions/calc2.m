function [runningFSWmedian, runningTSWmedian] = calc2(dataIn, thresholdIn, timespanIn )

timespan = timespanIn;  % 80 minutes
a.data = dataIn;
threshold = thresholdIn;

nRows = length(a.data);
startRow = 1;
endRow = nRows;

disp('START CALC2 -------------------------------------------------------')

runningFSWmedian = zeros( size( a.data(startRow:endRow) ));
runningFSWmedian(:,:) = NaN;

runningTSWmedian = zeros( size( a.data(startRow:endRow) ));
runningTSWmedian(:,:) = NaN;

% this is my index of TSW/FSW  TSW = 1/T; FSW = 0/F;
inTSW = zeros( size( a.data(startRow:endRow) ));
inTSW(:,:) = NaN;
inTSW = a.data >= threshold;
%added:
TSWdat = a.data(inTSW);
FSWdat = a.data(~inTSW);

%%
for iRow = startRow:endRow %1:57600;  %1:10000 %nRows

    % we have current position: iRow
    % x timestamps backwards: iStartBack
    % if we have passed the end of set timespan
    % i.e. not going to go out of range if we subtract timespan
    % from current location
    if iRow > timespan
        % we have passed timespan end point
        iStartBack = iRow - timespan;
    else
        % we are within first timespan, if we subtract we will go out of
        % range, so start with 1
        iStartBack = 1;
    end
 
    % subset data & index
     runData = a.data(iStartBack:iRow);
     runTSWIdx = inTSW(iStartBack:iRow);
    if (inTSW(iRow)) % if we are inTSW
        runningTSWmedian(iRow) = median(runData(runTSWIdx));
%         runningTSWmedian(iRow) = mean(TSWdat(iStartBack:iRow));
    else
        runningFSWmedian(iRow) = median(runData(~runTSWIdx));
%         runningFSWmedian(iRow) = mean(FSWdat(iStartBack:iRow));
    end
    
end   % end for loop of data records


disp('FINISH-------------------------------------------------------------')

end   % end function