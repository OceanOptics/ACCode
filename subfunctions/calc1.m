function [runningMinForw, runningMedForw] = calc1(dataIn, timespanIn )

timespan = timespanIn;  % 80 minutes
a.data = dataIn;

nRows = length(a.data);
startRow = 1;
endRow = nRows;

disp('START CALC 1 ------------------------------------------------------')

runningMinForw = zeros( size( a.data(startRow:endRow) ));
runningMinForw(:,:) = NaN;

runningMedForw = zeros( size( a.data(startRow:endRow) ));
runningMedForw(:,:) = NaN;

% create subset of a.data
%%
    for iRow = startRow:endRow 

        % set up data

        % we have current position: iRow
        % x timestamps backwards: iStartBack
        % x timestamps forward: iEndForward
        % if we have passed the end of set timespan
        % i.e. not going to go out of range if we subtract timespan
        % from current location

        if iRow > endRow - timespan
            % if we are within last timespan, if we add timespan, we will go
            % out of range
            iEndForward = endRow;
        else
            % we are able to add timespan to current row
            iEndForward = iRow + timespan;
        end

        runningdataforward = a.data(iRow:iEndForward);

        runningMinForw(iRow) = min(runningdataforward);
        runningMedForw(iRow) = median(runningdataforward);    



    end   % end for loop
disp('FINISH CALC 1 -----------------------------------------------------')
end   % end function