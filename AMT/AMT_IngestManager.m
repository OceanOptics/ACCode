%% IngestManager - SCRIPT which takes a set of declared variables and loads 
% all the necessary data for processing AC data for a single YEAR_DAY.
% For each of the file types, IngestManager calls the FileLoader for that
% type.  The FileLoader creates the appropriate data structures and returns
% them to the IngestManager.
% IngestManager only needs to know which file types are available for a
% given cruise.
%
% Requires:
%    ACDeviceFile - Description
%    FlowFileLoader - Description
%    TSGFileLoader - Description
%    ACFileLoader - opens data files using specified script and creates
%    data objects.
%
% Outputs:
%    allData - name of the data object this script creates and saves to
%    disk.  This is an "ACPlusAncillary" data object.
%    output2 - Description
%
%
% Other m-files required: none
% Subfunctions: readIngestParameters, log4m, getFilesNextPrev
% MAT-files required: none
%
% See also: ACDeviceFile, OceanLoggerFileLoader, SeatexFileLoader, 
% ACFileLoader, ACPlusAncillaryData

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% Sep 2016


%----------------------------- BEGIN CODE ---------------------------------

%% Create file names
% make a directory for this year day, if it doesn't already exist
if ~exist(fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY))
    mkdir(fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY))
end

matFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, params.INGEST.DATA_OUTPUT_FILE);
paramsFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, 'params');
%% Create device file object
L.debug('IngestManager', 'creating ACDeviceFile object');
devFile = ACDeviceFile( params.INGEST.DEVICE_FILE_LOCATION, params.INGEST.DEVICE_FILE_TYPE );

%% Start to ingest files
% For each of the main data file types for this cruise, get a list of files for
% the specified yearday, the last file from the previous day, and the
% first file from the next day.  
% Call the appropriate Loader for the file type (i.e. FlowFileLoader) and
% get it to create a data structure(s) for the appropriate data contained
% in each file type.


%%
% generate list of AC files
% needs to be customized for different file naming conventions:
% AC_FILE_FORMAT = 'amt24_14_%03d_24_ACS.*;

accurr = sprintf(params.INGEST.AC_FILE_FORMAT, params.INGEST.YEAR_DAY);
acnext = sprintf(params.INGEST.AC_FILE_FORMAT, params.INGEST.YEAR_DAY+1);
acprev = sprintf(params.INGEST.AC_FILE_FORMAT, params.INGEST.YEAR_DAY-1);
acFiles = getFilesNextPrev(params.INGEST.AC_DIRECTORY, accurr, acnext, acprev);

prepacsOut = fullfile(params.INGEST.PREPACS_OUTPUT_DIRECTORY);
%%
% load data from AC files
% ACFileLoader(fileNameListIn, importmethodIn, deviceFileIn, outputLocationIn)
acFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, params.INGEST.AC_ONLY_OUTPUT_FILE);
if params.INGEST.USE_ACFILELOADER == true
    if strcmp( params.INGEST.AC_FILE_LOADER_TYPE, 'bin')
        acfl = ACFileLoader( acFiles, params.INGEST.AC_IMPORT_METHOD_NAME, ...
            params.INGEST.DEVICE_FILE_LOCATION, prepacsOut, params.INGEST.AC_UNITS, ...
            devFile.NumberWavelengths, params.INGEST.PREPACS_BIN);
    elseif strcmp( params.INGEST.AC_FILE_LOADER_TYPE, 'dat')
        acfl = ACFileLoaderDat( acFiles, params.INGEST.AC_IMPORT_METHOD_NAME, ...
            params.INGEST.DEVICE_FILE_LOCATION, prepacsOut, params.INGEST.AC_UNITS, ...
            devFile.NumberWavelengths);
    else
        acfl = ACFileLoaderWAP( acFiles, params.INGEST.AC_IMPORT_METHOD_NAME, ...
            params.INGEST.DEVICE_FILE_LOCATION, prepacsOut, params.INGEST.AC_UNITS, ...
            devFile.NumberWavelengths);        
    end
    acData = acfl.loadData();
    
    % if we want to save this AC data to disk, to keep from having to
    % reprocess
    if params.INGEST.SAVE_AC_TO_DISK == true
        save( acFileName, 'acData');
    end
else   % not calling ACFileLoader
    acFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, params.INGEST.AC_ONLY_OUTPUT_FILE);
    acData = load(acFileName);
end
%% GET ANCILLARY DATA
% For AMT24, the OCEANLOGGER data files contain:
%    date
%    salinity
%    flowrate
%    temperature
% and the SEATEX files contain:
%    latitude & longitude

%% generate list of OceanLogger files
% needs to be customized for different file naming conventions:
ol_curr = sprintf(params.INGEST.OCEANLOGGER_FILE_FORMAT, params.INGEST.YEAR_DAY);
ol_next = sprintf(params.INGEST.OCEANLOGGER_FILE_FORMAT, params.INGEST.YEAR_DAY+1);
ol_prev = sprintf(params.INGEST.OCEANLOGGER_FILE_FORMAT, params.INGEST.YEAR_DAY-1);

% call script getFilesNextPrev, return a list of files
olFiles = getFilesNextPrev(params.INGEST.OCEANLOGGER_DIRECTORY, ol_curr, ol_next, ol_prev);

%% load data from OceanLogger File, using list of files generated above
% 
if size(olFiles) >= 1
    olfl = OceanLoggerFileLoader( olFiles, params.INGEST.OCEANLOGGER_IMPORT_METHOD_NAME, ...
        'temperature', params.INGEST.TEMP_UNITS, 'salinity', params.INGEST.SAL_UNITS, ...
        'flow', params.INGEST.FLOW_UNITS );
    oceanloggerData = olfl.loadData();
end;

%% generate list of SEATEX files
seatex_curr = sprintf(params.INGEST.SEATEX_FILE_FORMAT, params.INGEST.YEAR_DAY);
seatex_next = sprintf(params.INGEST.SEATEX_FILE_FORMAT, params.INGEST.YEAR_DAY+1);
seatex_prev = sprintf(params.INGEST.SEATEX_FILE_FORMAT, params.INGEST.YEAR_DAY-1);
% call script getFilesNextPrev, return a list of files
seatexFiles = getFilesNextPrev(params.INGEST.OCEANLOGGER_DIRECTORY, ...
    seatex_curr, seatex_next, seatex_prev);

%% load data from SEATEX File, using list of files generated above
% 
if size(seatexFiles) >= 1
    seaxtexfl = SEATEXFileLoader( seatexFiles, ...
        params.INGEST.SEATEX_IMPORT_METHOD_NAME, 'gps', params.INGEST.GPS_UNITS );
    seatexData = seaxtexfl.loadData();
end;
%%
if exist('acData','var')
    FileTypes = { oceanloggerData, acData, seatexData};
end;
nDataSources = length(FileTypes);

%% Plot
% loop through each of the data objects from each of the data files
% plot data from same datafiles together on same plot
for iDataSources = 1:nDataSources;
    
    DataTypes = FileTypes{1,iDataSources};
    nDataTypes = length(DataTypes);    

    for iDataTypes = 1:nDataTypes;
   
       tempDataType = DataTypes{iDataTypes};
       fprintf('Data Type Name: %s\n', tempDataType.Name)
       fprintf('Data Type: %s\n', tempDataType.Type)
   
       figure(iDataSources); 
       subplot(nDataTypes, 1, iDataTypes);
       tempDataType.plotData();
       fprintf('min: %s\n', datestr(min(tempDataType.DataObject.Time), 'dd-mm-yyyy HH:MM:SS'))
       fprintf('max: %s\n', datestr(max(tempDataType.DataObject.Time), 'dd-mm-yyyy HH:MM:SS'))
       hold on
       saveas(gcf,  strcat(matFileName, tempDataType.Type));
   end
end

%% Load data into ACPlusAncillaryData object
dataMatrix = [seatexData, oceanloggerData, acData];
allData = ACPlusAncillaryData( devFile, dataMatrix, params );

%% Save data to disk
save( matFileName, 'allData');
save( paramsFileName, 'params');

%% Make final plot
figure(4)
hold on;
grid on;
ax1 = subplot(3,1,1);

scatter(allData.aData.DataObject.Time, allData.aData.DataObject.Data(:,20), 'k');
dynamicDateTicks;
xlabel('Time');
ylabel('a Data');
legend('a Data');
if params.INGEST.VALVE_EXISTS
    ax2 = subplot(3,1,2);
    dynamicDateTicks;
    plot(allData.ValveData.DataObject.Time, allData.ValveData.DataObject.Data, 'c');
    xlabel('Time');
    ylabel('Valve Data');
    legend('Valve Data');
end;
ax3 = subplot(3,1,3);
dynamicDateTicks;
scatter(allData.TemperatureData.DataObject.Time, allData.TemperatureData.DataObject.Data(:,1));
xlabel('Time');
ylabel('Temperature Data');
legend('Temperature Data');
if params.INGEST.VALVE_EXISTS

    linkaxes([ax1, ax2, ax3], 'x');
else
        linkaxes([ax1, ax3], 'x');
end

%% close variables

if params.INGEST.CLEAR_VARS
    clear accurr;
    clear acData;
    clear acFileName;
    clear acFiles;
    clear acfl;
    clear acnext;
    clear acprev;
    clear ans;
    clear ax1;
    clear ax2;
    clear ax3;
    clear dataMatrix;
    clear DataTypes;
    clear devFile;
    clear ffl;
    clear FileTypes;
    clear flowcurr;
    clear flowFiles;
    clear flowData;
    clear flowfiles;
    clear flownext;
    clear flowprev;
    clear iDataSources;
    clear iDataTypes;
    clear matFileName;
    clear nDataSources;
    clear nDataTypes;
    clear paramsFileName;
    clear prepacsOut;
    clear tempDataType;
    clear tfl;
    clear tsgcurr;
    clear tsgData;
    clear tsgData;
    clear tsgFiles;
    clear tsgnext;
    clear tsgprev;
end;
%----------------------------- END CODE ---------------------------------