function printTransitionTimes(startTimesIn, endTimesIn, fileNameIn )
    % create a file of start and end times
    startTimes = startTimesIn;
    endTimes = endTimesIn;
    fileName = fileNameIn;
        
    fid = fopen(fileName, 'wt');
    [r,c] = size(startTimes);
    if r<c
        startTimes = startTimes';
    end;
    for i = 1:size(startTimes)
        fprintf(fid, 'start,%s, 1\n', datestr(startTimes(i)'));
    end;
    [r,c] = size(endTimes);
    if r<c
        endTimes = endTimes';
    end;
    for i = 1:size(endTimes)
        fprintf(fid, 'end,%s, 1\n', datestr(endTimes(i)'));
    end;
    fclose(fid);
    clear fileName;
end