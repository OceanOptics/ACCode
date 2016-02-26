%% ProcessingManager - SCRIPT which takes pre-processed data, saved as an 
% ACPlusAncillaryData object called "allData" and does the main data
% processing for a and c data, as well as binning temperature, salinity and
% GPS data (for Latitude and Longitude).
%
% STEP 0:   Create ProcessedData Object
% STEP 1A:  Process Filtered (FSW): (LEVEL L3 (BINNED) DATA)
% STEP 1B:  Process Total (TSW)
% STEP 2:   Get Filtered Spectra (INTERPOLATE): (LEVEL L4 (FILTERED) DATA)
% STEP 3:   Calculate Particulate: (LEVEL L5: PARTICULATE)
% STEP 4:   Compute new uncertainty for p
% STEP 5&6: Remove WL over 750 & interpolate a to c (Level L6: BELOW750)
% STEP 7:   Correct mismatch in spectral band positions between a and c 
%                                                   (Level L7: MATCHEDWL)
% STEP 8:   Scattering/Residual Temperature Correction(s) (LEVEL L8:
% CORRECTED)
% STEP 9:   Correct attenuation
% STEP 10:  Calculate final uncertainty for ap
% STEP 11:  Spectral unsmoothing: (LEVEL L9: UNSMOOTH)
% STEP 12:  Bin TSG data
%
% Requires: the .mat file saved in PreProcessingManager, if not in memory
%
% Outputs:  a ProcessedData object, also saved as a mat file
%
% Other m-files required: readIngestParameters, ACPlusAncillaryData,
% ProcessedData, ACData, ProcessedData, TemperatureData, SalinityData,
% GPSData
% Subfunctions: spectral_unsmooth, ResidTempScatCorr, AttTempCorr
% MAT-files required: "acsPREPROC_(YEAR)_(DAY)" saved by PreProcessingManager
%
% See also: readIngestParameters, ACPlusAncillaryData, ProcessedData,
% spectral_unsmooth, ResidTempScatCorr, AttTempCorr, ACData,
% TemperatureData, SalinityData, GPSData
% 
% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% May 2015; Last revision: 13-08-15
% 

%----------------------------- BEGIN CODE ---------------------------------
%%Load data file from disk
% Set to TRUE if opening raw data file from disk
% Set to FALSE if running automatically (from script) or you've been
% running it manually and have just run IngestManager and data objects are
% still in memory.

if params.RUN.LOAD_PREPROCESS_DATA_FROM_DISK
    matFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
        strcat('acsPREPROC', '_', num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY)));
    load(matFileName);
%     paramsFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'params');
%     load(paramsFileName);    
end;

% %%
% % LOGGING SETUP 
% % create log file 
% L = log4m.getLogger(params.LOG_FILE);
% 
% % set output level
% L.setCommandWindowLevel(L.DEBUG);
% L.setLogLevel(L.DEBUG);

%% Use the indices created by .separateFSWTSW to actually separate data
%cFSW
cFSW = zeros(size(allData.cData.SyncDataObject.Data));
cFSW(:) = NaN;
cFSW(allData.cData.FSWIndex,:) = allData.cData.SyncDataObject.Data(allData.cData.FSWIndex,:);

% cTSW
cTSW = zeros(size(allData.cData.SyncDataObject.Data));
cTSW(:) = NaN;
cTSW(allData.cData.TSWIndex,:) = allData.cData.SyncDataObject.Data(allData.cData.TSWIndex,:);

%aFSW
aFSW = zeros(size(allData.aData.SyncDataObject.Data));
aFSW(:) = NaN;
aFSW(allData.aData.FSWIndex,:) = allData.aData.SyncDataObject.Data(allData.aData.FSWIndex,:);

% aTSW
aTSW = zeros(size(allData.aData.SyncDataObject.Data));
aTSW(:) = NaN;
aTSW(allData.aData.TSWIndex,:) = allData.aData.SyncDataObject.Data(allData.aData.TSWIndex,:);

%% COMMENT THIS SECTION IN TO: make shorter data for testing
% cFSW = cFSW(1:50000,:);
% cTSW = cTSW(1:50000,:);
% aFSW = aFSW(1:50000,:);
% aTSW = aTSW(1:50000,:);
% 
% rawA = allData.aData.SyncDataObject.Data(1:50000,:);
% rawC = allData.cData.SyncDataObject.Data(1:50000,:);
% rawAT = allData.aData.SyncDataObject.Time(1:50000,:);
% rawCT = allData.cData.SyncDataObject.Time(1:50000,:);


% THIS ISN'T REALLY RAW THOUGH -- IT'S SYNCED...
rawA = allData.aData.SyncDataObject.Data(:,:);
rawC = allData.cData.SyncDataObject.Data(:,:);
rawAT = allData.aData.SyncDataObject.Time(:,:);
rawCT = allData.cData.SyncDataObject.Time(:,:);

%% STEP 0: create ProcessedData Object
pd = ProcessedData( allData.DeviceFile, params );

% -----------Need a lot of logic here re. a/c existing -------------------
% Create Level 1 Data ('raw')
%  function obj = setVar( obj, varNameIn, levelNameIn, dataNameIn, dataIn )
pd.setVar('a', 'raw', 'data', rawA); %allData.aData.SyncDataObject.Data );
pd.setVar('c', 'raw', 'data', rawC); %allData.cData.SyncDataObject.Data );
pd.setVar('a', 'raw', 'timestamps', rawAT); %allData.aData.SyncDataObject.Time );
pd.setVar('c', 'raw', 'timestamps', rawCT); %allData.cData.SyncDataObject.Time );

% Create Level 2 Data ('preprocessed')
pd.setVar('aTSW', 'preprocessed', 'data', aTSW);
pd.setVar('aTSW', 'preprocessed', 'timestamps', rawAT); %allData.aData.SyncDataObject.Time );
pd.setVar('aFSW', 'preprocessed', 'data', aFSW);
pd.setVar('aFSW', 'preprocessed', 'timestamps', rawAT); %allData.aData.SyncDataObject.Time );
pd.setVar('cTSW', 'preprocessed', 'data', cTSW);
pd.setVar('cTSW', 'preprocessed', 'timestamps', rawCT); %allData.cData.SyncDataObject.Time );
pd.setVar('cFSW', 'preprocessed', 'data', cFSW);
pd.setVar('cFSW', 'preprocessed', 'timestamps', rawCT); %allData.cData.SyncDataObject.Time );

pd.aExists = true;
pd.cExists = true;

%% ------------------------------------------------------------------------
% STEP 1A: PROCESS FILTERED DATA (FSW): CREATES LEVEL L3 (BINNED) DATA
% STEP 1A: PROCESS 'TOTAL' DATA (TSW)
% -------------------------------------------------------------------------
% 
% 1.  Bin it
% for each filtered period
% create minute bins
% find data between 2.5 and 97.5%
% flag data outside these as suspect
% calculate median and mean for rest
% if | median - mean | > (0.0005 + (0.05*median))
%   then flag
% rebin unflagged data to 1 minute
% calculate mean, median, std
% calculate percentiles

% CREATE MINUTE BINS (total and filtered)
pd.processBins();

% calculate suspect data
pd.calcSuspectData();
% pd.calcSuspectFSWData()
% pd.calcSuspectFSWData();

% PLOT
% NEED NOTE HERE ABOUT WHAT IS BEING PLOTTED
% THESE ARE DEBUG PLOTS
if params.RUN.CREATE_DEBUG_PLOTS
    pd.plotFSWSuspectData(31, 20);
    pd.plotTSWSuspectData(32, 20);
end;
% ------------------------------------------------------------------------ 
% STEP 2: Get Filtered Spectra (INTERPOLATE): CREATES LEVEL L4 (FILTERED)
% DATA
% -------------------------------------------------------------------------
%  Interpolate between dissolved spectra using linear interpolation
%  Calculate STD for interpolated bins by finding (RANGE of start and end
%  filter bins/2)

% STEP 2A: INTERPOLATE BETWEEN DISSOLVED USING LINEAR INTERPOLATION
%
% STEP 2B: CALCULATE UNCERTAINTY for INTERPOLATED BINS
%
% This is (MAX - MIN)/2
pd.findFSWBinMedians();
pd.interpolateFiltered();

%%  Plot Data
%THIS IS DEBUG PLOT
pd.plotACInterpolatedData(33, 20);
saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
    strcat(num2str(params.INGEST.YEAR_DAY), '_filt_spec_median')));

%% ------------------------------------------------------------------------
% STEP 3: CALCULATE PARTICULATE: CREATES LEVEL L5: PARTICULATE
% STEP 4: COMPUTE NEW UNCERTAINTY FOR P
%
% % Subtract filtered measurements (FROM TOM)
% bin_cp = bin_rawc - cgi;
% bin_ap_uncorr = bin_rawa - agi;

% have cTSW and cFSW
pd.calcParticulate();

% check ap & cp timestamps are the same, otherwise make the same
ap_size = size(pd.var.ap.L5.timestamps);
cp_size = size(pd.var.cp.L5.timestamps);
if ap_size(1) ~= cp_size(1)
    L.info('ProcessingManager','ap size and cp size diff')
    if abs(cp_size(1) - ap_size(1)) == 1
        % only one timestamp in differnece - difference of seconds
        if cp_size(1) > ap_size(1)
            %cp is bigger
            % if end-1 for cp is same as end for ap
            if sum(datevec(pd.var.ap.L5.timestamps(end)) == datevec(pd.var.cp.L5.timestamps(end-1))) == 6
                pd.var.cp.L5.timestamps(end,:) = [];
                pd.var.cp.L5.data(end,:) = [];
                pd.var.cp.L5.uncertainty(end,:) = [];
            else
                L.error('ProcessingManager','gap between a tiimestamps and c timestamps is larger than 1');
            end;
        elseif ap_size(1) > cp_size(1)
            %ap is bigger
            % if end-1 for ap is same as end for cp
            if sum(datevec(pd.var.ap.L5.timestamps(end-1)) == datevec(pd.var.cp.L5.timestamps(end))) == 6
                pd.var.ap.L5.timestamps(end) = [];
                pd.var.ap.L5.data(end) = [];
                pd.var.ap.L5.uncertainty(end) = [];
            else
                L.error('ProcessingManager','gap between a tiimestamps and c timestamps is larger than 1');
            end;
        else
            L.error('ProcessingManager','neither a nor c seems to be bigger');
        end;
    else
        L.error('ProcessingManager','gap between a tiimestamps and c timestamps is larger than 1');
    end;
end;
% check again
ap_size2 = size(pd.var.ap.L5.timestamps);
cp_size2 = size(pd.var.cp.L5.timestamps);
if ap_size2(1) ~= cp_size2(1)
    L.error('ProcessingManager', 'tried to fix size gap between a and c tiemstamps and failed')
else
    L.info('ProcessingManager', 'adjusted gap between and c timestamps by 1');
end;

% plot particulate debug plots
if params.RUN.CREATE_DEBUG_PLOTS
    pd.plotCpVsTime(34);
    pd.plotApVsTime(35);
end;


%% STEP 5&6: REMOVE WAVELENGTHS OVER 750 & INTERPOLATE A 750 wl of data

if pd.cExists == true
    pd.removeWLAfter750('c');
else
    if pd.aExists == true
        % WHY DID I COMMENT THIS OUT?
%       pd.removeWLAfter750('a');
    else
        L.error('ProcessingManager', 'aExists and cExists both false');
    end;
end;



%% STEP 7: Correct mismatch in spectral band positions between a and c 
% Only can match the wavelengths, if we have both a and c data
% a will be interpolated to c wavelengths

% also outputs ag and cg
if pd.cExists == true
    if pd.aExists == true
        L.debug('ProcessingManager', 'cExists == true; aExists == true');
        pd.correctSpectralBandMismatch('c');
    else
        L.debug('ProcessingManager', 'cExists == true; aExists == false');
    end;
else
    L.debug('ProcessingManager', 'cExists == false');
end;
    

%% Step 8: Scattering/Residual Temperature Correction

if pd.cExists == true
    if pd.aExists == true
        L.debug('ProcessingManager', 'cExists == true; aExists == true');
        
        % more than one correction is possible
        if params.PROCESS.SCATTERING_CORR_SLADE
            pd.scatteringCorr('SLADE');
        end;
        if params.PROCESS.SCATTERING_CORR_ROTTGERS
            pd.scatteringCorr('ROTTGERS');
        end;
        if params.PROCESS.SCATTERING_CORR_FLAT
            pd.scatteringCorr('FLAT')
        end;
    else
        L.debug('ProcessingManager', 'cExists == true; aExists == false');
        % do nothing - go to correct attenuation
    end;
else
    L.debug('ProcessingManager', 'cExists == false');
    if pd.aExists == true
        L.debug('ProcessingManager', 'cExists == false; aExists == true');        
            pd.scatteringCorr('FLAT');
    end;
end;    
        
pd.plotApCorr(37, 'slade');
pd.plotApCorr(38, 'rottgers');

%% Step 9: CORRECT ATTENUATION
if pd.cExists == true
    if pd.aExists == true
        L.debug('ProcessingManager', 'cExists == true; aExists == true');
        pd.attenuationCorr('WITHA');
        
        if params.RUN.CREATE_DEBUG_PLOTS
            tsToUse = find( ~isnan(pd.var.cp.L8.data(:,1)), 1, 'first');

            figure(39)
            hold on; 
            grid on;
            plot( pd.var.cp.L8.wavelengths, pd.var.cp.L8.data(tsToUse,:), '*b')
            plot( pd.var.cp.L7.wavelengths, pd.var.cp.L7.data(tsToUse,:), '*c')
            legend('corrected cp', 'uncorrected cp')
            title('correcting c using a')
        end;
         
    else
        
        % DO NOT THINK THIS IS WORKING
        
        L.debug('ProcessingManager', 'cExists == true; aExists == false');
        pd.attenuationCorr('WITHOUTA');
        
        tsToUse = find( ~isnan(pd.var.cp.L9.data(:,1)), 1, 'first');
        if params.RUN.CREATE_DEBUG_PLOTS
            figure(39)
            hold on; 
            grid on;
            plot( pd.var.cp.L9.wavelengths, pd.var.cp.L9.data(tsToUse,:), '*b')
            plot( pd.var.cp.L8.wavelengths, pd.var.cp.L8.data(tsToUse,:), '*c')
            legend('corrected cp', 'uncorrected cp')
            title('correcting c without using a')

            figure(310)
            hold on; 
            grid on;
            plot( pd.var.cp.L9.wavelengths, pd.var.cp.L9.data(tsToUse,:), 'b')
            plot( pd.var.cp.L8.wavelengths, pd.var.cp.L8.data(tsToUse,:), 'c')
            legend('corrected cp', 'uncorrected cp')
            title('correcting c without using a')
        end;
    end;
else
    L.debug('ProcessingManager', 'no attenuation correction - no c data');
end;
% pd.attenuationCorr('WITHA');
%%  STEP 10: CALCULATE FINAL UNCERTAINTY FOR ap

% for ap:
% uncertainty in processing

pd.computeAPUncertaintyBetweenCorrections();

%% ------------------------------------------------------------------------
% STEP 11: SPECTRAL UNSMOOTHING: CREATES LEVEL L6: UNSMOOTH
%

if params.PROCESS.UNSMOOTH_DATA
    pd.unsmooth();
    
    if params.RUN.CREATE_DEBUG_PLOTS

        % DEBUG PLOT TO CHECK
        atsToUse = find( ~isnan(pd.var.ap.L8.data_slade(:,1)), 1, 'first');
        ctsToUse = find( ~isnan(pd.var.cp.L9.data(:,1)), 1, 'first');
        pd.plotUnsmoothVsSmooth(311, atsToUse, ctsToUse);
        
    end;
else
    L.info('ProcessingManager', 'Not doing unsmoothing');
end;


%% STEP 12: Bin TSG data

tData = allData.TemperatureData.DataObject.Data(:,1);
tTime = allData.TemperatureData.DataObject.Time;

% THIS HAS TO BE CHANGED:
% WHAT WAS THIS NOTE FOR??  BECAUSE i WAS CALLING THE DATA EXPLICITLY?
if pd.cExists == true
    acTime = pd.var.cp.L5.timestamps;
elseif pd.aExists == true
    acTime = pd.var.ap.L5.timestamps;
else
    L.error('ProcessingManager', 'neither cExists or aExists set to true');
end;

% bin data -- method is stored in params
allData.TemperatureData.bin(acTime, params.PROCESS);
allData.SalinityData.bin(acTime, params.PROCESS);
allData.GPSData.bin(acTime, params.PROCESS);

%% put all on one big plot
% THIS IS A GOOD INTERMEDIATE PLOT TO OUTPUT AND SAVE
% Temperature Data
figure(312)
ax1 = subplot(2,2,1);
grid on;
hold on;
scatter(allData.TemperatureData.var.L3.BinnedTimestamps, allData.TemperatureData.var.L3.BinnedLabTempData, 'b')
dynamicDateTicks;
xlabel('Binned Timestamps');
ylabel('Binned Temperature - C');
legend('Binned Instrument Temperature Data')

% Salinity Data
ax2 = subplot(2,2,2);
grid on;
hold on;
dynamicDateTicks;
plot(allData.SalinityData.var.L3.BinnedTimestamps, allData.SalinityData.var.L3.BinnedData, 'b')
xlabel('Binned Timestamps');
ylabel('Binned Salinity');
legend('Binned Salinity Data');


% GPS Data
ax3 = subplot(2,2,3);
grid on;
hold on;
scatter(allData.GPSData.var.L3.BinnedLatData, allData.GPSData.var.L3.BinnedLonData, [], allData.GPSData.var.L3.BinnedTimestamps)
xlabel('Latitude')
ylabel('Longitude')
legend('Binned GPS Data');

% ap Data
ax4 = subplot(2,2,4);
hold on;
grid on;
dynamicDateTicks;
plot(pd.var.cp.L5.timestamps, pd.var.cp.L8.data(:,60))
xlabel('Binned Timestamps');
ylabel('cp - corrected, unsmoothed data');

linkaxes([ax1,ax2,ax4], 'x');
saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, strcat(num2str(params.INGEST.YEAR_DAY), '_temp_sal_gps_ap')));

%% Save file
% WANTED TO SAVE PD, BUT ALSO, RESAVE ALLDATA??
    matFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
        strcat('acsPROC', '_', num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY)));

    save( '-v7.3', matFileName, 'pd');
% 
    matFileName2 = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
        strcat('acsPREPROC', '_', num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY)));

    save(matFileName2, 'allData');
%% close variables

if params.INGEST.CLEAR_VARS
    clear acTime;
    clear aFSW;
    clear ans;
    clear ap_size;
    clear ap_size2;
    clear atsToUse;
    clear aTSW;
    clear ax1;
    clear ax2;
    clear ax3;
    clear ax4;
    clear cFSW;
    clear cTSW;
    clear cp_size;
    clear cp_size2;
    clear ctsToUse;
    clear matFileName;
    clear matFileName2;
    clear rawA;
    clear rawAT;
    clear rawC;
    clear rawCT;
    clear tData;
    clear tsToUse;
    clear tTime;
end;