function [runningFSWmedian, runningTSWmedian] = findFSWTSWRunMeds( dataIn, timeIn, timespanIn, samplingFreqIn )

    % check these variables
    timespan = timespanIn; %19200 %80 minutes 2400 %timespanIn; %2400;(2400 readings == 10 min
    samplingFreq = samplingFreqIn; %4: 4*60 = 240  FOR MINUTE?

    disp('starting first pass')
  
    [runningMinForw, runningMedForw] = calc1(dataIn, timespan);   
    
    manualAdjustmentMins = timespan/samplingFreq; 
    manualAdjustmentSecs = manualAdjustmentMins*60;

%     figure(100);
%     hold on;
%     grid on;
%     plot(timeIn, dataIn);
%     dynamicDateTicks;
%     plot(timeIn + datenum(0,0,0,0,manualAdjustmentMins,0), runningMedForw, 'g*');
%     plot(timeIn + datenum(0,0,0,0,manualAdjustmentMins,0), runningMedForw - mean(runningMedForw - runningMinForw)/2, 'r*')
    % end plots
    
    % figure out new threshold:
   
    runningThresholdTmp = runningMedForw - mean(runningMedForw - runningMinForw)/2;
    runningThreshold = zeros(size(runningMedForw));
    runningThreshold(:,:) = nanmean(runningThresholdTmp(1:manualAdjustmentSecs));
    runningThreshold(manualAdjustmentSecs+1:end) = runningThresholdTmp(1:end-manualAdjustmentSecs);
    
    disp('starting second pass')
    [runningFSWmedian, runningTSWmedian] = calc2(dataIn, runningThreshold, timespan);
    
%     figure(101)
%     hold on
%     grid on
%     plot(timeIn, dataIn)
%     plot(timeIn, runningFSWmedian, 'r*')
%     plot(timeIn, runningTSWmedian, 'g*')
%     legend('data', 'running FSW median', 'running TSW median');
%     dynamicDateTicks
end