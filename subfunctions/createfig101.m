    figure(101)
    hold on
    grid on
    plot(allData.aData.DataObject.Time, allData.aData.DataObject.Data(:,20))
    plot(timeIn, allData.aData.runningFSWmedian, 'r*')
    plot(timeIn, allData.aData.runningTSWmedian, 'g*')
    legend('data', 'running FSW median', 'running TSW median');
    dynamicDateTicks