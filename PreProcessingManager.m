%% PreProcessingManager - SCRIPT which takes raw data, saved as an 
% ACPlusAncillaryData object called "allData" and gets it ready for 
% processing by the ProcessingManager.
% The PreProcessingManager has the following three main jobs:
% -- 1: Find good TSW/FSW periods in:
%           -- a/c data
%           -- flow/valve data
%       Methodology:
%           -- use flow and valve state data
%           -- use flow data only
%           -- use a/c data only
% -- 2: Synchronize:  a/c data with flow/valve data
% -- 3: Separate TSW/FSW data
%
% It also allows for a set offset to be applied to the start of each FSW
%
% Requires: the .mat file saved in IngestManager, if not in memory
%
% Outputs:  saves allData as a .mat file
%
% Other m-files required: readIngestParameters, ACPlusAncillaryData,
% Subfunctions: importTransitionFile, checkTransitionTimes
% MAT-files required: none
%
% See also: readIngestParameters, ACPlusAncillaryData,
% importTransitionFile, checkTransitionTimes

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% May 2015; Last revision: 13-08-15
% July/August 2016 - major revision in way transitions are managed

%--------------------------------- BEGIN CODE -----------------------------
%% Load data file from disk
matFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, params.INGEST.DATA_OUTPUT_FILE);
paramsFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'params');

if params.RUN.LOAD_INGEST_DATA_FROM_DISK
    L.info('PreProcessingManager', 'Loading ingest data from disk');
    load(matFileName);
else
    L.info('PreProcessingManager', 'Using ingest data in memory');
end;
clear matFileName;


%% Find intervals
if params.RUN.FIND_INITIAL_INTERVALS

    %% make shorter for testing:
    % a.data = a.data(1:10000,:);
    % a.time = a.time(1:10000,:);
    % c.data = c.data(1:10000,:);
    % c.time = c.time(1:10000,:);

    % Make an initial plot of a data, alongside valve data, if available
    % Set ylim manually if needed
    if params.RUN.CREATE_DEBUG_PLOTS
        figure(20)
        hold on
        grid on
        plot(allData.aData.DataObject.Time(:,:), allData.aData.DataObject.Data(:, 20), 'b');
        if params.INGEST.VALVE_EXISTS
            plot(allData.ValveData.DataObject.Time, allData.ValveData.DataObject.Data, 'g');
            legend('a data', 'valve data');
            title('a data & valve data - unsynchronized');
            % Set ylim here:
            %     ylim([-.20, .25])
        else
            legend('a data');
            title('a data - No valve data exists');
        end;
        % ylim([0,.2])
        dynamicDateTicks;
    end  %#params.RUN.CREATE_DEBUG_PLOTS

    %% if there is valve on/off information available, use it to mark the 
    % transition periods, rather than identifying them from the actual data
    if  params.PREPROCESS.USE_VALVE_DATA

        allData.ValveData.findTransitions();
        if params.RUN.CREATE_DEBUG_PLOTS
            allData.ValveData.plotTransitions();
        end;

        if params.PREPROCESS.USE_VALVE_FOR_FLOW_TRANSITIONS

            % set transition points for flow by using those from valve
            allData.FlowData.setTransitions( allData.ValveData.TransStartData, ...
                allData.ValveData.TransEndData, allData.ValveData.TransStartTime, ...
                allData.ValveData.TransEndTime )
            if params.RUN.CREATE_DEBUG_PLOTS
                allData.FlowData.plotTransitions()
            end;
        else
            L.info('PreProcessingManager', 'Not using valve to set flow transitions');
        end;  %if params.PREPROCESS.USE_VALVE_FOR_FLOW_TRANSITIONS
        
        if params.PREPROCESS.USE_VALVE_FOR_AC_TRANSITIONS
            % valve start/end times are for FSW
            one_min = datenum(0,0,0,0,1,0);
            two_min = datenum(0,0,0,0,2,0);
            thirty_secs = datenum(0,0,0,0,0,30);
            
            ValveOn = allData.ValveData.TransStartTime;
            ValveOff = allData.ValveData.TransEndTime;
            
            allData.aData.setGoodTransitions(ValveOn + two_min, ...
                ValveOff - thirty_secs, 'FSW', 'raw');
            allData.cData.setGoodTransitions(ValveOn + two_min, ...
                ValveOff - thirty_secs, 'FSW', 'raw');

            allData.aData.setGoodTransitions(ValveOff + two_min, ...
                ValveOn - one_min, 'TSW', 'raw');
            allData.cData.setGoodTransitions(ValveOff + two_min, ...
                ValveOn - one_min, 'TSW', 'raw');              
        else
            L.info('PreProcessingManager', 'Not using valve to set ac');
        end;
    else   % params.PREPROCESS.USE_VALVE_DATA
        L.info('PreProcessingManager', 'Valve not available.');
    end;
    if params.INGEST.FLOW_EXISTS
            L.info('PreProcessingManager', 'Flow exists');
        if params.PREPROCESS.USE_FLOW_DATA

            % smooth the data
            allData.FlowData.setSmoothData();
            if params.RUN.CREATE_DEBUG_PLOTS
                allData.FlowData.plotSmoothData();
            end;

            allData.FlowData.setRunningMeds(params.PREPROCESS.FLOW_PP_TIMESPAN, params.PREPROCESS.SAMPLING_FREQ)
            if params.RUN.CREATE_DEBUG_PLOTS
                figure(12)

                allData.FlowData.plotRunningMeds()
            end;

            % find transition times from running TSW median
            allData.FlowData.findTransitions()
            if params.RUN.CREATE_DEBUG_PLOTS
                figure(23)
                allData.FlowData.plotTransitions()
            end;
        else
            L.info('PreProcessingManager', 'Not using flow to find flow transitions');
        end;
        if params.PREPROCESS.USE_FLOW_FOR_AC_TRANSITIONS
            L.info('PreProcessingManager', 'Using flow to find ac transitions');
            
        else
            L.info('PreProcessingManager', 'Not using flow to find ac transitions');
        end;
    end;   %FLOW_EXISTS
    
    if params.PREPROCESS.USE_AC_FOR_AC_TRANSITIONS
            
        % -------------------------------------------------------------
        % a AND c DATA SECTION
        % -------------------------------------------------------------
        % -------------------------------------------------------------
        % Smooth a and c data, find the running medians and then use the running 
        % medians to find the transitions between FSW and TSW

        % a data

        allData.aData.setSmoothData();
        allData.aData.setRunningMeds(params.PREPROCESS.PP_TIMESPAN, ...
            params.PREPROCESS.SAMPLING_FREQ);   

        allData.aData.findAndSetTransitionPoints('TSW')
        allData.aData.findAndSetTransitionPoints('FSW')

        % plot TSW/FSW transitions to check
        if params.RUN.CREATE_DEBUG_PLOTS
            figure(21)
            hold on;
            grid on;
            allData.aData.plotSmoothData();
            allData.aData.plotRunningMeds();
            allData.aData.plotInitFSWTransitions();
            allData.aData.plotInitTSWTransitions();
        end;

        % c data

        allData.cData.setSmoothData();
        allData.cData.setRunningMeds(params.PREPROCESS.PP_TIMESPAN, ...
            params.PREPROCESS.SAMPLING_FREQ);
        allData.cData.findAndSetTransitionPoints('TSW')
        allData.cData.findAndSetTransitionPoints('FSW')

        % plot TSW/FSW transitions to check
        if params.RUN.CREATE_DEBUG_PLOTS
            figure(22)
            allData.cData.plotSmoothData();
            allData.cData.plotRunningMeds();
            allData.cData.plotInitFSWTransitions()
            allData.cData.plotInitTSWTransitions()
        end;

    end;  % IF USE_FLOW_FOR_AC_TRANSITIONS
    
    % save current as mat file
    matFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
        strcat('acsPREPROC_PARTIAL', '_', num2str(params.INGEST.YEAR), ...
        '_', num2str(params.INGEST.YEAR_DAY)));
    save( matFileName, 'allData');
    clear matFileName;

else  % if params.RUN.FIND_INITIAL_INTERVALS
    
    L.info('PreProcessingManager', 'Skipping finding intitial intervals');
    matFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
        strcat('acsPREPROC_PARTIAL', '_', num2str(params.INGEST.YEAR), ...
        '_', num2str(params.INGEST.YEAR_DAY)));
    load(matFileName);   
    clear matFileName;
end % if FIND_INITIAL_INTERVALS


%% CHECKING TRANSITIONS 
% find first good transition for flow, a/c
if params.PREPROCESS.USE_FLOW_DATA || params.PREPROCESS.USE_VALVE_FOR_FLOW_TRANSITIONS
    L.info('PreProcessingManager','Calling allData.FlowData.checkTransitions()');
    %FSWFreq, FSWDuration
    allData.FlowData.checkTransitions(params.PREPROCESS.TSW_DURATION, ...
        params.PREPROCESS.FSW_DURATION)
else
    L.info('PreProcessingManager','No flow data to check transitions on');
end;

%% If manually editing, but using checked/filtered transitions, go ahead
%  and check/filter them
% If running manually, and want to check/filter transitions first
% OR not running manually
%
if (params.RUN.MANUAL_MODE && params.RUN.MANUALLY_EDIT_CHECKED_TRANSITIONS) ...
        || (~params.RUN.MANUAL_MODE && ~params.RUN.REPROCESS_MANUAL)
    
    Transitions = allData.aData.getTransitions('raw');

    % check and set aFSW transitions ----------------------------------------------
    [ goodStartsOut, goodEndsOut ] = ...
        checkTransitionTimes( Transitions.FSWtransitions.StartTime,...
        Transitions.FSWtransitions.EndTime, allData.aData.DataObject.Time );
                            
    [ goodStartTimes, goodEndTimes ] = filterTransitionTimes( goodStartsOut, ...
        goodEndsOut, params.PREPROCESS.FSW_DURATION, ...
        params.PREPROCESS.FSW_DUR_TOLERANCE, params.PREPROCESS.CYCLE_FREQ, ...
        params.PREPROCESS.FREQ_TOLERANCE );

    allData.aData.setGoodTransitions(goodStartTimes, goodEndTimes, 'FSW','preprocessed');
    
    % check and set aTSW transitions ----------------------------------------------
        [ goodStartsOut, goodEndsOut ] = ...
        checkTransitionTimes( Transitions.TSWtransitions.StartTime,...
        Transitions.TSWtransitions.EndTime, allData.aData.DataObject.Time );
                            
    [ goodStartTimes, goodEndTimes ] = filterTransitionTimes( goodStartsOut, ...
        goodEndsOut, params.PREPROCESS.TSW_DURATION, ...
        params.PREPROCESS.TSW_DUR_TOLERANCE, params.PREPROCESS.CYCLE_FREQ, ...
        params.PREPROCESS.FREQ_TOLERANCE );

    allData.aData.setGoodTransitions(goodStartTimes, goodEndTimes, 'TSW', 'preprocessed');
  
    % check cFSW transitions -----------------------------------------------
    
    Transitions = allData.cData.getTransitions('raw');
    % check and set cFSW transitions ----------------------------------------------
    [ goodStartsOut, goodEndsOut ] = ...
        checkTransitionTimes( Transitions.FSWtransitions.StartTime,...
        Transitions.FSWtransitions.EndTime, allData.aData.DataObject.Time );
                            
    [ goodStartTimes, goodEndTimes ] = filterTransitionTimes( goodStartsOut, ...
        goodEndsOut, params.PREPROCESS.FSW_DURATION, ...
        params.PREPROCESS.FSW_DUR_TOLERANCE, params.PREPROCESS.CYCLE_FREQ, ...
        params.PREPROCESS.FREQ_TOLERANCE );

    allData.cData.setGoodTransitions(goodStartTimes, goodEndTimes, 'FSW', 'preprocessed');
    
    % check and set cTSW transitions ----------------------------------------------
        [ goodStartsOut, goodEndsOut ] = ...
        checkTransitionTimes( Transitions.TSWtransitions.StartTime,...
        Transitions.TSWtransitions.EndTime, allData.aData.DataObject.Time );
                            
    [ goodStartTimes, goodEndTimes ] = filterTransitionTimes( goodStartsOut, ...
        goodEndsOut, params.PREPROCESS.TSW_DURATION, ...
        params.PREPROCESS.TSW_DUR_TOLERANCE, params.PREPROCESS.CYCLE_FREQ, ...
        params.PREPROCESS.FREQ_TOLERANCE );

    allData.cData.setGoodTransitions(goodStartTimes, goodEndTimes, 'TSW', 'preprocessed');
    
else
    L.info('PreProcessingManager','Not using automatic transtions, using manually-set');
end %~MANUAL_MODE


%% stop here for manual processing, if we're editing files
if params.RUN.MANUAL_MODE && ~params.RUN.REPROCESS_MANUAL
    L.info('PreProcessingManager', 'manual mode');
    return;
else
    L.info('PreProcessingManager', 'not in manual mode');
end;

%% ------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%  START MANUAL PROCESSING HERE:
%  1. CREATE PLOT TO CHECK c TRANSITIONS - plot raw and possibly edited
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
if params.INGEST.VALVE_EXISTS
    figure(24)
    hold on;
    grid on;
    allData.cData.plotData();
    plot(allData.ValveData.DataObject.Time, allData.ValveData.DataObject.Data, 'c');
    if params.RUN.MANUALLY_EDIT_CHECKED_TRANSITIONS
        allData.cData.plotTSWGoodTransitions();
        legend('c data', 'valve data', 'identified checked TSW periods - start', ...
        'identified checked TSW periods - end');
    else
        allData.cData.plotInitTSWTransitions();
        legend('c data', 'valve data', 'identified initial TSW periods - start', ...
        'identified initial TSW periods - end');        
    end;
else
    figure(24)
    hold on;
    grid on;
    allData.cData.plotData();
    if params.RUN.MANUALLY_EDIT_CHECKED_TRANSITIONS
        allData.cData.plotTSWGoodTransitions();
        legend('c data', 'identified checked TSW periods - start', ...
        'identified checked TSW periods - end');
    else
        allData.cData.plotInitTSWTransitions();
        legend('c data', 'identified initial TSW periods - start', ...
        'identified initial TSW periods - end');        
    end;
end;
saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, strcat(num2str(params.INGEST.YEAR_DAY), '_cTSW_transitions')));

%% <-------------MANUAL PROCESSING CONTINUES HERE------------------------->
if params.RUN.MANUAL_MODE || params.RUN.REPROCESS_MANUAL
    %
    % DON'T RUN FROM HERE, RUN FROM SECTION BELOW
    %
    %% -------- SECTION 1 - c TSW - PART 1A-----------------
    % 2.  "RUN AND ADVANCE" FROM HERE
    % -----------------------------------------------------
    
    % if we're not reprocessing this - run this first part:
    if ~params.RUN.REPROCESS_MANUAL

        % PROCESS C TSWs
        % create a file of start and end times
        if params.RUN.MANUALLY_EDIT_RAW_TRANSITIONS
            Transitions = allData.cData.getTransitions('raw');
        elseif params.RUN.MANUALLY_EDIT_CHECKED_TRANSITIONS
            Transitions = allData.cData.getTransitions('preprocessed');
        else
            L.error('PreProcessingManager', 'Need type of transitions');
        end;

        fileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'transitions_cTSW.txt');
        printTransitionTimes(Transitions.TSWtransitions.StartTime, ...
            Transitions.TSWtransitions.EndTime, fileName);
        
        L.info('PreProcessingManager', 'Reprocessing MANUAL_MODE, cTSW file ready');
    else
        L.info('PreProcessingManager', 'Reprocessing MANUAL_MODE, not recreating transitions text file');
    end %~REPROCESS_MANUAL
    
    %% ------- SECTION1 - c TSW - PART 1B -----------
    % 3.  "RUN AND ADVANCE" FROM HERE AFTER YOU HAVE:
    % -----------------------------------------------
    % A.  In your PROCESSED directory, open up transitions_cTSW.txt
    % B.  Delete any unwanted transitions, using Figure 20 as a guide
    % 4.  Save the file, in the same directory as transtions_cTSW2.txt
    % 5.  "RUN AND ADVANCE" 
        
    %% ------- SECTION1 - c TSW - PART 2-----------------
    % 4.  "RUN AND ADVANCE" FROM HERE
    % ---------------------------------------------------
    
    % open file back up
    filename = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'transitions_cTSW2.txt');
    [type, timestamp, flag] = importTransitionFile(filename);

    % manipulate data out of file
    startIndex = strcmpi(type,'start');
    cTSWTransitions.startTimestamps = timestamp(startIndex);
    cTSWTransitions.startFlags = flag(startIndex);
    cTSWTransitions.endTimestamps = timestamp(~startIndex);
    cTSWTransitions.endFlags = flag(~startIndex);

    % check transitions
    [ goodStarts, goodEnds ] = ...
        checkTransitionTimes(cTSWTransitions.startTimestamps,...
        cTSWTransitions.endTimestamps, allData.aData.DataObject.Time );

    % NOT FILTERING BECAUSE SET MANUALLY

    % call method to assign good starts
    allData.cData.setGoodTransitions(goodStarts, goodEnds, 'TSW', 'preprocessed');

    L.info('PreProcessingManager', 'Reprocessing MANUAL_MODE, done with cTSW');
    %%
    if params.INGEST.FLOW_EXISTS
        figure(25)
        hold on;
        grid on;
        allData.cData.plotData();
        plot(allData.ValveData.DataObject.Time, allData.ValveData.DataObject.Data, 'c');
        if params.RUN.MANUALLY_EDIT_CHECKED_TRANSITIONS
            allData.cData.plotFSWGoodTransitions();
            legend('c data', 'valve data', 'identified checked FSW periods - start', ...
            'identified checked FSW periods - end');
        else
            allData.cData.plotInitFSWTransitions();
            legend('c data', 'valve data', 'identified initial FSW periods - start', ...
            'identified initial FSW periods - end');        
        end;
    else
        figure(25)
        hold on;
        grid on;
        allData.cData.plotData();
        if params.RUN.MANUALLY_EDIT_CHECKED_TRANSITIONS
            allData.cData.plotFSWGoodTransitions();
            legend('c data', 'identified checked FSW periods - start', ...
            'identified checked FSW periods - end');
        else
            allData.cData.plotInitFSWTransitions();
            legend('c data', 'identified initial FSW periods - start', ...
            'identified initial FSW periods - end');        
        end;
    end;
    saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, strcat(num2str(params.INGEST.YEAR_DAY), '_cFSW_transitions')));

    %% ------- SECTION 2 - c FSW- PART 1A-----------------
    % 5.  "RUN AND ADVANCE" FROM HERE
    % ---------------------------------------------------
    if ~params.RUN.REPROCESS_MANUAL
    
        % PROCESS c FSWs
        if params.RUN.MANUALLY_EDIT_RAW_TRANSITIONS
            Transitions = allData.cData.getTransitions('raw');
        elseif params.RUN.MANUALLY_EDIT_CHECKED_TRANSITIONS
            Transitions = allData.cData.getTransitions('preprocessed');
        else
            L.error('PreProcessingManager', 'Need type of transitions');
        end;        
        % create a file of start and end times
        fileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'transitions_cFSW.txt');
        printTransitionTimes(Transitions.FSWtransitions.StartTime, ...
            Transitions.FSWtransitions.EndTime, fileName);
        L.info('PreProcessingManager', 'Reprocessing MANUAL_MODE, cFSW file ready');
    else
        L.info('PreProcessingManager', 'Reprocessing MANUAL_MODE, not recreating transitions text file');
    end; %~REPROCESS_MANUAL
    
    %% ------- SECTION 2 - c FSW - PART 1B -----------
    % 6.  "RUN AND ADVANCE" FROM HERE AFTER YOU HAVE:
    % ------------------------------------------------
        % A.  In your PROCESSED directory, open up transitions_cFSW.txt
        % B.  Delete any unwanted transitions, using Figure 20 as a guide
        % C.  Save the file, in the same directory as transtions_cFSW2.txt
        % D.  "RUN AND ADVANCE" 
        
    %% ------- SECTION 2 - c FSW - PART 2-----------------
    % 7.  "RUN AND ADVANCE" FROM HERE
    % ----------------------------------------------------
    
    % open file back up
    filename = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'transitions_cFSW2.txt');
    [type, timestamp, flag] = importTransitionFile(filename);

    % manipulate data out of file
    startIndex = strcmpi(type,'start');
    cFSWTransitions.startTimestamps = timestamp(startIndex);
    cFSWTransitions.startFlags = flag(startIndex);
    cFSWTransitions.endTimestamps = timestamp(~startIndex);
    cFSWTransitions.endFlags = flag(~startIndex);

    % check transitions
    [ goodStarts, goodEnds ] = ...
        checkTransitionTimes( cFSWTransitions.startTimestamps,...
        cFSWTransitions.endTimestamps, allData.aData.DataObject.Time );

    % call method to assign good ends
    allData.cData.setGoodTransitions(goodStarts, goodEnds, 'FSW', 'preprocessed');

    L.info('PreProcessingManager', 'Reprocessing MANUAL_MODE, done with cFSW');

    %% ---------------------------------------------------------------
    % ------------------------ PROCESS A -----------------------------
    %  8. CREATE PLOT TO CHECK a TRANSITIONS
    if params.INGEST.FLOW_EXISTS
        figure(25)
        hold on;
        grid on;
        allData.aData.plotData();
        plot(allData.ValveData.DataObject.Time, allData.ValveData.DataObject.Data, 'c');
        allData.aData.plotInitFSWTransitions();
        allData.aData.plotInitTSWTransitions();
        if params.RUN.MANUALLY_EDIT_CHECKED_TRANSITIONS
            allData.aData.plotFSWGoodTransitions();
            allData.aData.plotTSWGoodTransitions();
        end;
        legend('a data', 'valve data', 'identified FSW periods - start', ...
            'identified FSW periods - end',  ...
            'identified TSW periods - start', ...
            'identified TSW periods - end',...
            'good FSW transition start', 'good FSW transition end', ...
            'good TSW transition start', 'good TSW transition end');
    else
        figure(25)
        allData.aData.plotData();
        allData.aData.plotInitFSWTransitions();
        allData.aData.plotInitTSWTransitions();
        if params.RUN.MANUALLY_EDIT_CHECKED_TRANSITIONS
            allData.aData.plotFSWGoodTransitions();
            allData.aData.plotTSWGoodTransitions();
        end;
        legend('a data', 'identified FSW periods - start', ...
            'identified FSW periods - end',  ...
            'identified TSW periods - start', ...
            'identified TSW periods - end',...
            'good FSW transition start', 'good FSW transition end', ...
            'good TSW transition start', 'good TSW transition end');
    end;    

    
    %% --------SECTION 3 - a TSW - PART 1A-----------------
    % 9.  "RUN AND ADVANCE" FROM HERE
    % ---------------------------------------------------
        
     if ~params.RUN.REPROCESS_MANUAL 
         if params.RUN.MANUALLY_EDIT_RAW_TRANSITIONS
            Transitions = allData.aData.getTransitions('raw');
        elseif params.RUN.MANUALLY_EDIT_CHECKED_TRANSITIONS
            Transitions = allData.aData.getTransitions('preprocessed');
        else
            L.error('PreProcessingManager', 'Need type of transitions');
        end;
         
         
        % create a file of start and end times
        fileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'transitions_aTSW.txt');
        printTransitionTimes(Transitions.TSWtransitions.StartTime, ...
            Transitions.TSWtransitions.EndTime, fileName);
        L.info('PreProcessingManager', 'Reprocessing MANUAL_MODE, aTSW file ready');

     else
             
        L.info('PreProcessingManager', 'Reprocessing MANUAL_MODE, not recreating transitions text file');
     end;
    %% ------- SECTION 3 - a TSW - PART 1B -----------
    % 10.  "RUN AND ADVANCE" FROM HERE AFTER YOU HAVE:
    % -----------------------------------------------
        % A.  In your PROCESSED directory, open up transitions_aTSW.txt
        % B.  Delete any unwanted transitions, using Figure 20 as a guide
        % 4.  Save the file, in the same directory as transtions_aTSW2.txt
        % 5.  "RUN AND ADVANCE" 
        
    %% ------- SECTION 3 - a TSW - PART 2-----------------
    % 11.  "RUN AND ADVANCE" FROM HERE
    % ---------------------------------------------------
    filename = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'transitions_aTSW2.txt');
    [type, timestamp, flag] = importTransitionFile(filename);

    % manipulate data out of file
    startIndex = strcmpi(type,'start');
    aTSWTransitions.startTimestamps = timestamp(startIndex);
    aTSWTransitions.startFlags = flag(startIndex);
    aTSWTransitions.endTimestamps = timestamp(~startIndex);
    aTSWTransitions.endFlags = flag(~startIndex);

    % check transitions
    [ goodStarts, goodEnds ] = ...
        checkTransitionTimes( aTSWTransitions.startTimestamps,...
        aTSWTransitions.endTimestamps, allData.aData.DataObject.Time );

    % call method to assign good starts
    allData.aData.setGoodTransitions(goodStarts, goodEnds, 'TSW', 'preprocessed');

    %% ------- SECTION 4 - a FSW- PART 1A-----------------
    % 12.  "RUN AND ADVANCE" FROM HERE
    % ---------------------------------------------------
    if ~params.RUN.REPROCESS_MANUAL
         if params.RUN.MANUALLY_EDIT_RAW_TRANSITIONS
            Transitions = allData.aData.getTransitions('raw');
        elseif params.RUN.MANUALLY_EDIT_CHECKED_TRANSITIONS
            Transitions = allData.aData.getTransitions('preprocessed');
        else
            L.error('PreProcessingManager', 'Need type of transitions');
        end;
        % create a file of start and end times
        fileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'transitions_aFSW.txt');
        printTransitionTimes(Transitions.FSWtransitions.StartTime, ...
            Transitions.FSWtransitions.EndTime, fileName);
        L.info('PreProcessingManager', 'Reprocessing MANUAL_MODE, aFSW file ready');

    else
            
        L.info('PreProcessingManager', 'Reprocessing MANUAL_MODE, not recreating transitions text file');
    end;
    %% ------- SECTION 4 - a FSW - PART 1B -----------
    % 13.  "RUN AND ADVANCE" FROM HERE AFTER YOU HAVE:
    % ------------------------------------------------
        % A.  In your PROCESSED directory, open up transitions_aFSW.txt
        % B.  Delete any unwanted transitions, using Figure 20 as a guide
        % C.  Save the file, in the same directory as transtions_aFSW2.txt
        % D.  "RUN AND ADVANCE" 
        
    %% ------- SECTION 4 - a FSW - PART 2-----------------
    % 13.  "RUN AND ADVANCE" FROM HERE
    % ----------------------------------------------------
    filename = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'transitions_aFSW2.txt');
    [type, timestamp, flag] = importTransitionFile(filename);

    % manipulate data out of file
    startIndex = strcmpi(type,'start');
    aFSWTransitions.startTimestamps = timestamp(startIndex);
    aFSWTransitions.startFlags = flag(startIndex);
    aFSWTransitions.endTimestamps = timestamp(~startIndex);
    aFSWTransitions.endFlags = flag(~startIndex);

    % check transitions
    [ goodStarts, goodEnds ] = ...
        checkTransitionTimes( aFSWTransitions.startTimestamps,...
        aFSWTransitions.endTimestamps, allData.aData.DataObject.Time );

    % call method to assign good starts & ends
    allData.aData.setGoodTransitions(goodStarts, goodEnds, 'FSW', 'preprocessed');

end;  %if params.RUN.MANUAL_MODE

% <-----------------MANUAL PROCESSING ENDS HERE--------------------------->
%  IF YOU HAVE BEEN RUNNING MANUALLY, RUN AND ADVANCE THROUGH EACH SECTION
%  BELOW:
% <----------------------------------------------------------------------->

%%
% plot to show how a and c data are lined up
if params.RUN.CREATE_DEBUG_PLOTS
    if params.INGEST.VALVE_EXISTS
        figure(26)
        hold on;
        grid on;
        allData.aData.plotData();
        plot(allData.ValveData.DataObject.Time, allData.ValveData.DataObject.Data, 'c');
        allData.aData.plotInitFSWTransitions();
        allData.aData.plotInitTSWTransitions();
        allData.aData.plotFSWGoodTransitions();
        allData.aData.plotTSWGoodTransitions();
        legend('a data', 'valve data', 'identified FSW periods - start', ...
            'identified FSW periods - end',  ...
            'identified TSW periods - start', ...
            'identified TSW periods - end',...
            'good FSW transition start', 'good FSW transition end', ...
            'good TSW transition start', 'good TSW transition end');
        else
        figure(26)
        allData.aData.plotData();
        allData.aData.plotInitFSWTransitions();
        allData.aData.plotInitTSWTransitions();
        allData.aData.plotFSWGoodTransitions();
        allData.aData.plotTSWGoodTransitions();
        legend('a data', 'identified FSW periods - start', ...
            'identified FSW periods - end',  ...
            'identified TSW periods - start', ...
            'identified TSW periods - end',...
            'good FSW transition start', 'good FSW transition end', ...
            'good TSW transition start', 'good TSW transition end');
    end;
end;  %params.RUN.CREATE_DEBUG_PLOTS

if params.RUN.CREATE_DEBUG_PLOTS
    if params.INGEST.VALVE_EXISTS
        figure(27)
        hold on;
        grid on;
        allData.cData.plotData();
        plot(allData.ValveData.DataObject.Time, allData.ValveData.DataObject.Data, 'c');
        allData.cData.plotInitFSWTransitions();
        allData.cData.plotInitTSWTransitions();
        allData.cData.plotFSWGoodTransitions();
        allData.cData.plotTSWGoodTransitions();
        legend('c data', 'valve data', 'identified FSW periods - start', ...
            'identified FSW periods - end',  ...
            'identified TSW periods - start', ...
            'identified TSW periods - end',...
            'good FSW transition start', 'good FSW transition end', ...
            'good TSW transition start', 'good TSW transition end');
    else
        figure(27)
        allData.cData.plotData();
        allData.cData.plotInitFSWTransitions();
        allData.cData.plotInitTSWTransitions();
        allData.cData.plotFSWGoodTransitions();
        allData.cData.plotTSWGoodTransitions();
        legend('c data', 'identified FSW periods - start', ...
            'identified FSW periods - end',  ...
            'identified TSW periods - start', ...
            'identified TSW periods - end',...
            'good FSW transition start', 'good FSW transition end', ...
            'good TSW transition start', 'good TSW transition end');
    end;
end; %if params.RUN.CREATE_DEBUG_PLOTS

%% SYNCHRONIZE Data

if params.INGEST.FLOW_EXISTS && params.PREPROCESS.SYNC_AC_TO_FLOW
L.info('PreProcessingManager', 'Flow data exists for syncing to');
    allData.aData.syncTo(allData.FlowData);
    allData.cData.syncTo(allData.FlowData);
else
    L.info('PreProcessingManager', 'Either No flow exists for sync, OR we aren''t using it.  Copying regular data into SyncDataObject');
    allData.aData.SyncDataObject = allData.aData.DataObject;
    allData.cData.SyncDataObject = allData.cData.DataObject;
end;

%% apply offset set in parameters

if params.PREPROCESS.OFFSET_FSW_STARTS
    L.info('PreProcessingManager', 'Applying offset to FSW start as set in params');
    allData.aData.offsetFSWStartTimes(params.PREPROCESS.OFFSET_FSW_TIME);
    allData.cData.offsetFSWStartTimes(params.PREPROCESS.OFFSET_FSW_TIME);
else
    L.info('PreProcessingManager', 'Not applying FSW start offset');
end;

%% plot to check synchronization, and filter/total starts and ends
% THIS IS MAJOR QA/QC PLOT TO KEEP
figure(28)
title( sprintf('%s %s', params.INGEST.CRUISE, params.INGEST.CRUISE_LEG ))

ax1 = subplot(2,1,1);
plot(allData.aData.SyncDataObject.Time, allData.aData.SyncDataObject.Data(:,20), 'k');
hold on;
grid on;
dynamicDateTicks;
allData.aData.plotFSWGoodTransitions();
allData.aData.plotTSWGoodTransitions();
if params.INGEST.VALVE_EXISTS
    plot(allData.ValveData.DataObject.Time, allData.ValveData.DataObject.Data, 'c');
%     ylim([0,.15])
end
legend('sync a data', 'good fsw start', 'good fsw end' , 'good tsw start' , 'good tsw end', 'valve state')

ax2 = subplot(2,1,2);
plot(allData.cData.SyncDataObject.Time, allData.cData.SyncDataObject.Data(:,20), 'k');
hold on;
grid on;
dynamicDateTicks;
allData.cData.plotFSWGoodTransitions();
allData.cData.plotTSWGoodTransitions();
if params.INGEST.VALVE_EXISTS
    plot(allData.ValveData.DataObject.Time, allData.ValveData.DataObject.Data, 'c');
%     ylim([0,.6])
end
if params.INGEST.VALVE_EXISTS
    legend('sync c data', 'good fsw start', 'good fsw end' , 'good tsw start' , 'good tsw end', 'valve state')
else
    legend('sync c data', 'good fsw start', 'good fsw end' , 'good tsw start' , 'good tsw end')
end;
linkaxes([ax1, ax2], 'x');
saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, strcat(num2str(params.INGEST.YEAR_DAY), '_final_transitions')));

%% --------------------------------------
% separate TSW from FSW
% want two filters (TSW, FSW) and flag the transition data as suspect

allData.cData.separateTSWFSW();

%% plot TO CHECK DATA IS BEING SEPARATED CORRECTLY
% THIS IS A DEBUG PLOT
if params.RUN.CREATE_DEBUG_PLOTS
    figure(29);
    ax1 = subplot(3,1,1);
    hold on;
    grid on;
    allData.cData.plotSyncData()
    allData.cData.plotFSWGoodTransitions()
    allData.cData.plotTSWGoodTransitions()
    legend('c data', 'start FSW', 'end FSW', 'start TSW', 'end TSW')
    title( sprintf('%s %s %4u-%02u-%2u (yearday %3u)', ...
        params.INGEST.CRUISE, params.INGEST.CRUISE_LEG, params.INGEST.YEAR, ...
        params.INGEST.MONTH, params.INGEST.DAY, params.INGEST.YEAR_DAY))
    dynamicDateTicks;

    ax2 = subplot(3,1,2);
    hold on;
    grid on;
    plot(allData.cData.SyncDataObject.Time(allData.cData.FSWIndex), allData.cData.SyncDataObject.Data(allData.cData.FSWIndex), 'g');
    legend('Filtered c data')
    dynamicDateTicks;

    ax3 = subplot(3,1,3);
    hold on;
    grid on;
    plot(allData.cData.SyncDataObject.Time(allData.cData.TSWIndex), allData.cData.SyncDataObject.Data(allData.cData.TSWIndex), 'b');
    legend('Total c data')
    dynamicDateTicks;

    linkaxes([ax1, ax2, ax3], 'x');
    saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, strcat(num2str(params.INGEST.YEAR_DAY), '_fig29')));
end;  %if params.RUN.CREATE_DEBUG_PLOTS
%% - DO NOW FOR A
allData.aData.separateTSWFSW();

%% DEBUG PLOT
% plot TO CHECK DATA IS BEING SEPARATED CORRECTLY FOR A
% THIS IS A DEBUG PLOT
if params.RUN.CREATE_DEBUG_PLOTS
    figure(210);
    ax1 = subplot(3,1,1);
    hold on;
    grid on;
    allData.aData.plotSyncData()
    allData.aData.plotFSWGoodTransitions()
    allData.aData.plotTSWGoodTransitions()
    legend('a data', 'start FSW', 'end FSW', 'start TSW', 'end TSW')
    title( sprintf('%s %s %4u-%02u-%2u (yearday %3u)', ...
        params.INGEST.CRUISE, params.INGEST.CRUISE_LEG, params.INGEST.YEAR, params.INGEST.MONTH, params.INGEST.DAY, params.INGEST.YEAR_DAY))
    dynamicDateTicks;

    ax2 = subplot(3,1,2);
    hold on;
    grid on;
    plot(allData.aData.SyncDataObject.Time(allData.aData.FSWIndex), allData.aData.SyncDataObject.Data(allData.aData.FSWIndex), 'g');
    legend('Filtered a data')
    dynamicDateTicks;

    ax3 = subplot(3,1,3);
    hold on;
    grid on;
    plot(allData.aData.SyncDataObject.Time(allData.aData.TSWIndex), allData.aData.SyncDataObject.Data(allData.aData.TSWIndex), 'b');
    legend('Total a data')
    dynamicDateTicks;

    linkaxes([ax1,ax2,ax3], 'x');
    saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, strcat(num2str(params.INGEST.YEAR_DAY), '_fig210')));
end %if params.RUN.CREATE_DEBUG_PLOTS
%%
% save current as mat file
if params.RUN.SAVE_DATA
    matFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
        strcat('acsPREPROC', '_', num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY)));

    save( matFileName, 'allData');
    save( paramsFileName, 'params');
    clear matFileName;
    clear paramsFileName;
end;
%% close variables

if params.INGEST.CLEAR_VARS
    clear aFSWTransitions;
    clear aTSWTransitions;
    clear ans;
    clear ax1;
    clear ax2;
    clear goodEndsOut;
    clear goodEndTimes;
    clear goodStartsOut;
    clear goodStartTimes;
    clear goodEnds;
    clear goodStarts;
    clear cTSWTransitions;
    clear cFSWTransitions;
    clear filename;
    clear fileName;
    clear flag;
    clear startIndex;
    clear timestamp;
    clear Transitions;
    clear type;
end;

L.info('PreProcessingManager', 'FINISHED');

%--------------------------------- END CODE -----------------------------