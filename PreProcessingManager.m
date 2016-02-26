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

%--------------------------------- BEGIN CODE -----------------------------

matFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, params.INGEST.DATA_OUTPUT_FILE);
paramsFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'params');

if params.RUN.LOAD_INGEST_DATA_FROM_DISK
    L.info('PreProcessingManager', 'Loading ingest data from disk');
    load(matFileName);
else
    L.info('PreProcessingManager', 'Using ingest data in memory');
end;


%% Load data file from disk
% Set to TRUE if opening raw data file from disk
% Set to FALSE if running automatically (from script) or you've been
% running it manually and have just run IngestManager and data objects are
% still in memory.
if params.RUN.FIND_INITIAL_INTERVALS

    %% make shorter for testing:
    % a.data = a.data(1:10000,:);
    % a.time = a.time(1:10000,:);
    % c.data = c.data(1:10000,:);
    % c.time = c.time(1:10000,:);
    %% ----------------------------------------------------------------------
    % a AND c DATA SECTION
    % -----------------------------------------------------------------------
    % Make an initial plot of a data, alongside valve data, if available
    % Set ylim manually if needed
    if params.RUN.CREATE_DEBUG_PLOTS
        figure(20)
        hold on
        grid on
        plot(allData.aData.DataObject.Time(:,:), allData.aData.DataObject.Data(:, 20), 'b');
        if params.INGEST.FLOW_EXISTS
            plot(allData.ValveData.DataObject.Time, allData.ValveData.DataObject.Data, 'g');
            legend('a data', 'valve data');
            title('a data & valve data - unsynchronized');
        %     ylim([-.20, .25])
        else
            legend('a data');
            title('a data - No valve data exists');
        end;
        % ylim([0,.2])
        dynamicDateTicks;
    end  %#params.RUN.CREATE_DEBUG_PLOTS

    %% ----------------------------------------------------------------------
    % Smooth a and c data, find the running medians and then ue the running 
    % medians to find the transitions between FSW and TSW

    allData.aData.setSmoothData();
    allData.aData.setRunningMeds();   

    % find FSW Transitions
    allData.aData.findTransitions();

    % find TSW Transitions
    allData.aData.findTSWTransitions();

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

    %% c data
    % smooth data, find the running medians and then find the transitions
    % between FSW and TSW
    allData.cData.setSmoothData();
    allData.cData.setRunningMeds();
    allData.cData.findTransitions()
    allData.cData.findTSWTransitions();

    % plot TSW/FSW transitions to check
    if params.RUN.CREATE_DEBUG_PLOTS
        figure(22)
        allData.cData.plotSmoothData();
        allData.cData.plotRunningMeds();
        allData.cData.plotInitFSWTransitions()
        allData.cData.plotInitTSWTransitions()
    end;
    %% ----------------------------------------------------------------------
    % FLOW AND VALVE SECTION
    % -----------------------------------------------------------------------

    %% if there is valve on/off information available, use it to mark the 
    % transition periods, rather than identifying them from the actual data

    if  params.PREPROCESS.VALVE_AVAILABLE

        allData.ValveData.findTransitions();
        if params.RUN.CREATE_DEBUG_PLOTS
            allData.ValveData.plotTransitions();
        end;

        if params.PREPROCESS.USE_VALVE_AND_FLOW

            % set transition points for flow by using those from valve
            allData.FlowData.setTransitions( allData.ValveData.TransStartData, ...
                allData.ValveData.TransEndData, allData.ValveData.TransStartTime, allData.ValveData.TransEndTime )
            if params.RUN.CREATE_DEBUG_PLOTS
                allData.FlowData.plotTransitions()
            end;

        else
            L.info('PreProcessingManager', 'Not using valve or flow data');
        end   % end if USE_VALVE_AND_FLOW
    else
        L.info('PreProcessingManager','No valve data available');
    end   % end if VALVE_AVAILABLE

    %% ----------------------------------------------------------------------- 
    % find transitions in FLOW -- WITHOUT valve data
    % ------------------------------------------------------------------------

    if params.PREPROCESS.USE_FLOW_ONLY

        % smooth the data
        allData.FlowData.setSmoothData();
        if params.RUN.CREATE_DEBUG_PLOTS
            allData.FlowData.plotSmoothData();
        end;

        allData.FlowData.setRunningMeds()
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

    end   % if USE_FLOW_ONLY

    % save current as mat file
    matFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
        strcat('acsPREPROC_PARTIAL', '_', num2str(params.INGEST.YEAR), ...
        '_', num2str(params.INGEST.YEAR_DAY)));
    save( matFileName, 'allData');
    
else
    
    L.info('PreProcessingManager', 'Skipping finding intitial intervals');
    matFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
        strcat('acsPREPROC_PARTIAL', '_', num2str(params.INGEST.YEAR), ...
        '_', num2str(params.INGEST.YEAR_DAY)));
    load(matFileName);   
end % if FIND_INITIAL_INTERVALS

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
%  1. CREATE PLOT TO CHECK c TRANSITIONS
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
if params.INGEST.FLOW_EXISTS
    figure(24)
    hold on;
    grid on;
    allData.cData.plotData();
    plot(allData.ValveData.DataObject.Time, allData.ValveData.DataObject.Data, 'c');
    allData.cData.plotInitFSWTransitions();
    allData.cData.plotInitTSWTransitions();
    legend('c data', 'valve data','identified FSW periods - start', ...
        'identified FSW periods - end',  ...
        'identified TSW periods - start', ...
        'identified TSW periods - end');
%     linkaxes(p23,'x');
    saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, strcat(num2str(params.INGEST.YEAR_DAY), '_c_transitions')));
else
    figure(24)
    allData.cData.plotData();
    allData.cData.plotInitFSWTransitions();
    allData.cData.plotInitTSWTransitions();
    legend('c data', 'identified FSW periods - start', ...
        'identified FSW periods - end',  ...
        'identified TSW periods - start', ...
        'identified TSW periods - end');
    saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, strcat(num2str(params.INGEST.YEAR_DAY), '_c_transitions')));
end;

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
        
        fileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'transitions_cTSW.txt');
        fid = fopen(fileName, 'wt');
        [r,c] = size(allData.cData.TSWStartTime);
        if r<c
            allData.cData.TSWStartTime = allData.cData.TSWStartTime';
        end;
        for i = 1:size(allData.cData.TSWStartTime)
            fprintf(fid, 'start,%s, 1\n', datestr(allData.cData.TSWStartTime(i)'));
        end;
        [r,c] = size(allData.cData.TSWEndTime);
        if r<c
            allData.cData.TSWEndTime = allData.cData.TSWEndTime';
        end;
        for i = 1:size(allData.cData.TSWEndTime)
            fprintf(fid, 'end,%s, 1\n', datestr(allData.cData.TSWEndTime(i)'));
        end;
        fclose(fid);
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
    [ checkedStartTimes, checkedEndTimes, startFlagsOut, endFlagsOut ] = ...
        checkTransitionTimes( cTSWTransitions.startTimestamps, cTSWTransitions.endTimestamps, ...
        cTSWTransitions.startFlags,  cTSWTransitions.endFlags, ...
        allData.cData.DataObject.Time ) ;

    goodStartTimes = checkedStartTimes(startFlagsOut == 1);
    goodEndTimes = checkedEndTimes(endFlagsOut == 1);

    % call method to assign good starts
    allData.cData.setGoodTSWTransitions(goodStartTimes, goodEndTimes);

    L.info('PreProcessingManager', 'Reprocessing MANUAL_MODE, done with cTSW');
        
    %% ------- SECTION 2 - c FSW- PART 1A-----------------
    % 5.  "RUN AND ADVANCE" FROM HERE
    % ---------------------------------------------------
    if ~params.RUN.REPROCESS_MANUAL
    
        % PROCESS c FSWs
        
        % create a file of start and end times
        fileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'transitions_cFSW.txt');
        fid = fopen(fileName, 'wt');
        
        [r,c] = size(allData.cData.TransStartTime);
        if r<c
            allData.cData.TransStartTime = allData.cData.TransStartTime';
        end;
        for i = 1:size(allData.cData.TransStartTime)
            fprintf(fid, 'start,%s, 1\n', datestr(allData.cData.TransStartTime(i)));
        end;
        
        [r,c] = size(allData.cData.TransEndTime);
        if r<c
            allData.cData.TransEndTime = allData.cData.TransEndTime';
        end;        
        for i = 1:size(allData.cData.TransEndTime)
            fprintf(fid, 'end,%s, 1\n', datestr(allData.cData.TransEndTime(i)));
        end;
        fclose(fid);
        L.info('PreProcessingManager', 'Reprocessing MANUAL_MODE, created cFSW file');
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
    [ checkedStartTimes, checkedEndTimes, startFlagsOut, endFlagsOut ] = ...
        checkTransitionTimes( cFSWTransitions.startTimestamps, cFSWTransitions.endTimestamps, ...
        cFSWTransitions.startFlags,  cFSWTransitions.endFlags, ...
        allData.cData.DataObject.Time ); 

    goodStartTimes = checkedStartTimes(startFlagsOut == 1);
    goodEndTimes = checkedEndTimes(endFlagsOut == 1);

    % call method to assign good ends
    allData.cData.setGoodFSWTransitions(goodStartTimes, goodEndTimes);

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
        legend('a data', 'valve data','identified FSW periods - start', ...
            'identified FSW periods - end',  ...
            'identified TSW periods - start', ...
            'identified TSW periods - end');
    %     linkaxes(p23,'x');
    else
        figure(25)
        allData.aData.plotData();
        allData.aData.plotInitFSWTransitions();
        allData.aData.plotInitTSWTransitions();
        legend('a data', 'identified FSW periods - start', ...
            'identified FSW periods - end',  ...
            'identified TSW periods - start', ...
            'identified TSW periods - end');
    end;
    
    %% --------SECTION 3 - a TSW - PART 1A-----------------
    % 9.  "RUN AND ADVANCE" FROM HERE
    % ---------------------------------------------------
        
     if ~params.RUN.REPROCESS_MANUAL   
        % create a file of start and end times
        fileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'transitions_aTSW.txt');
        fid = fopen(fileName, 'wt');
        
        [r,c] = size(allData.aData.TSWStartTime);
        if r<c
            allData.aData.TSWStartTime = allData.aData.TSWStartTime';
        end;
        
        for i = 1:size(allData.aData.TSWStartTime)  %i = 1:size(allData.aData.TSWStartTime')
            fprintf(fid, 'start,%s, 1\n', datestr(allData.aData.TSWStartTime(i))');
        end;
        
        [r,c] = size(allData.aData.TSWEndTime);
        if r<c
            allData.aData.TSWEndTime = allData.aData.TSWEndTime';
        end;
        
        for i = 1:size(allData.aData.TSWEndTime)
            fprintf(fid, 'end,%s, 1\n', datestr(allData.aData.TSWEndTime(i))');
        end;
        fclose(fid);
        L.info('PreProcessingManager', 'Reprocessing MANUAL_MODE, creating transitions text file');

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
    [ checkedStartTimes, checkedEndTimes, startFlagsOut, endFlagsOut ] = ...
        checkTransitionTimes( aTSWTransitions.startTimestamps, aTSWTransitions.endTimestamps, ...
        aTSWTransitions.startFlags,  aTSWTransitions.endFlags, ...
        allData.aData.DataObject.Time ) ;

    goodStartTimes = checkedStartTimes(startFlagsOut == 1);
    goodEndTimes = checkedEndTimes(endFlagsOut == 1);

    % call method to assign good starts
    allData.aData.setGoodTSWTransitions(goodStartTimes, goodEndTimes);

    %% ------- SECTION 4 - a FSW- PART 1A-----------------
    % 12.  "RUN AND ADVANCE" FROM HERE
    % ---------------------------------------------------
    if ~params.RUN.REPROCESS_MANUAL
        
        % create a file of start and end times
        fileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'transitions_aFSW.txt');
        fid = fopen(fileName, 'wt');
        
        [r,c] = size(allData.aData.TransStartTime);
        if r<c
            allData.aData.TransStartTime = allData.aData.TransStartTime';
        end;
        
        for i = 1:size(allData.aData.TransStartTime)
            fprintf(fid, 'start,%s, 1\n', datestr(allData.aData.TransStartTime(i)));
        end;
        
        [r,c] = size(allData.aData.TransEndTime);
        if r<c
            allData.aData.TransEndTime = allData.aData.TransEndTime';
        end;
        
        for i = 1:size(allData.aData.TransEndTime)
            fprintf(fid, 'end,%s, 1\n', datestr(allData.aData.TransEndTime(i)));
        end;
        fclose(fid);
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
    [ checkedStartTimes, checkedEndTimes, startFlagsOut, endFlagsOut ] = ...
        checkTransitionTimes( aFSWTransitions.startTimestamps, aFSWTransitions.endTimestamps, ...
        aFSWTransitions.startFlags,  aFSWTransitions.endFlags, ...
        allData.aData.DataObject.Time ) ;

    goodStartTimes = checkedStartTimes(startFlagsOut == 1);
    goodEndTimes = checkedEndTimes(endFlagsOut == 1);

    % call method to assign good starts & ends
    allData.aData.setGoodFSWTransitions(goodStartTimes, goodEndTimes);        

    % WHAT ABOUT FLAGS? 

end;  %if params.RUN.MANUAL_MODE

% <-----------------MANUAL PROCESSING ENDS HERE--------------------------->
%  IF YOU HAVE BEEN RUNNING MANUALLY, RUN AND ADVANCE THROUGH EACH SECTION
%  BELOW:
% <----------------------------------------------------------------------->
%%
% find first good transition for flow, a/c
if params.INGEST.FLOW_EXISTS
    L.info('PreProcessingManager','Calling allData.FlowData.checkTransitions()');
    allData.FlowData.checkTransitions()
else
    L.info('PreProcessingManager','No flow data to check transitions on');
end;

% DON'T CALL THESE IF DONE MANUALLY or IF REPROCESSING MANUAL
if ~params.RUN.MANUAL_MODE && ~params.RUN.REPROCESS_MANUAL
%     L.info('PreProcessingManager','Calling allData.cData.checkTSWTransitions()');
%     allData.cData.checkTSWTransitions()
% 
%     L.info('PreProcessingManager','Calling allData.aData.checkTSWTransitions()');
%     allData.aData.checkTSWTransitions()
% 
%     allData.cData.checkFSWTransitions()
%     allData.aData.checkFSWTransitions()

    % check first a transitions - FSW
    [ goodStartsOut, goodEndsOut, startFlagsOut, endFlagsOut ] = ...
        checkTransitionTimes( allData.aData.TransStartTime, ... %rawStartsIn, 
                                allData.aData.TransEndTime, ... %rawEndsIn, 
                                allData.aData.TransStartFlag, ...%startFlagsIn, 
                                allData.aData.TransEndFlag, ...%endFlagsIn, 
                                allData.aData.DataObject.Time );
                            
    goodStartTimes = goodStartsOut(startFlagsOut == 1);
    goodEndTimes = goodEndsOut(endFlagsOut == 1);
    
    % call method to assign good starts
    allData.aData.setGoodFSWTransitions(goodStartTimes, goodEndTimes);
    
    % check a - TSW
    [ goodStartsOut, goodEndsOut, startFlagsOut, endFlagsOut ] = ...
        checkTransitionTimes( allData.aData.TSWStartTime, ... %rawStartsIn, 
                                allData.aData.TSWEndTime, ... %rawEndsIn, 
                                allData.aData.TSWStartFlag, ...%startFlagsIn, 
                                allData.aData.TSWEndFlag, ...%endFlagsIn, 
                                allData.aData.DataObject.Time ); 
                            
    goodStartTimes = goodStartsOut(startFlagsOut == 1);
    goodEndTimes = goodEndsOut(endFlagsOut == 1);
    
    
    % call method to assign good starts
    allData.aData.setGoodTSWTransitions(goodStartTimes, goodEndTimes);
    
        % next check c transitions - FSW
    [ goodStartsOut, goodEndsOut, startFlagsOut, endFlagsOut ] = ...
        checkTransitionTimes( allData.cData.TransStartTime, ... %rawStartsIn, 
                                allData.cData.TransEndTime, ... %rawEndsIn, 
                                allData.cData.TransStartFlag, ...%startFlagsIn, 
                                allData.cData.TransEndFlag, ...%endFlagsIn, 
                                allData.cData.DataObject.Time ) ;
                            
    goodStartTimes = goodStartsOut(startFlagsOut == 1);
    goodEndTimes = goodEndsOut(endFlagsOut == 1);

    % call method to assign good starts
    allData.cData.setGoodFSWTransitions(goodStartTimes, goodEndTimes);
    
           % next check c transitions - TSW
    [ goodStartsOut, goodEndsOut, startFlagsOut, endFlagsOut ] = ...
        checkTransitionTimes( allData.cData.TSWStartTime, ... %rawStartsIn, 
                                allData.cData.TSWEndTime, ... %rawEndsIn, 
                                allData.cData.TSWStartFlag, ...%startFlagsIn, 
                                allData.cData.TSWEndFlag, ...%endFlagsIn, 
                                allData.cData.DataObject.Time ); 
                            
    goodStartTimes = goodStartsOut(startFlagsOut == 1);
    goodEndTimes = goodEndsOut(endFlagsOut == 1);

    % call method to assign good starts
    allData.cData.setGoodTSWTransitions(goodStartTimes, goodEndTimes); 
else
    L.info('PreProcessingManager','Not using automatic transtions, using manually-set');
end %~MANUAL_MODE

%%
% plot to show how a and c data are lined up
if params.RUN.CREATE_DEBUG_PLOTS
    if params.INGEST.FLOW_EXISTS
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
    if params.INGEST.FLOW_EXISTS
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

if params.INGEST.FLOW_EXISTS
L.info('PreProcessingManager', 'Flow data exists for syncing to');
    allData.aData.syncTo(allData.FlowData);
    allData.cData.syncTo(allData.FlowData);
else
    L.info('PreProcessingManager', 'No flow exists for sync.  Copying regular data into SyncDataObject');
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
if params.INGEST.FLOW_EXISTS
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
if params.INGEST.FLOW_EXISTS
    plot(allData.ValveData.DataObject.Time, allData.ValveData.DataObject.Data, 'c');
%     ylim([0,.6])
end
if params.INGEST.FLOW_EXISTS
    legend('sync c data', 'good fsw start', 'good fsw end' , 'good tsw start' , 'good tsw end', 'valve state')
else
    legend('sync c data', 'good fsw start', 'good fsw end' , 'good tsw start' , 'good tsw end')
end;
linkaxes([ax1, ax2], 'x');
saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, strcat(num2str(params.INGEST.YEAR_DAY), '_fig28')));

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
matFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
    strcat('acsPREPROC', '_', num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY)));

save( matFileName, 'allData');
save( paramsFileName, 'params');

%% close variables

if params.INGEST.CLEAR_VARS
    clear ans;
    clear ax1;
    clear ax2;
    clear endFlagsOut;
%     clear f;
%     clear flow;
    clear goodEndsOut;
    clear goodEndTimes;
    clear goodStartsOut;
    clear goodStartTimes;
    clear matFileName;
    clear paramsFileName;
    clear startFlagsOut;
%     clear valve;
%     clear valveData;
%     clear valveTime;
    clear aFSWTransitions;
    clear aTSWTransitions;
    clear cFSWTransitions;
    clear fid
    clear filename;
    clear fileName;
    clear flag;
    clear startIndex;
    clear r;
    clear timestamp;
    clear type;
end;

L.info('PreProcessingManager', 'FINISHED');

%--------------------------------- END CODE -----------------------------