%% testNAAMESTSGFileLoader
% Script to test the NAAMES_TSGFileLoader object

disp('Starting TestNAAMESTSGFileLoader')

close all
clear all
%% ------------------------------------------------------------------------
% Get variable names
params = readIngestParameters();

%% ------------------------------------------------------------------------
% Query the directory for the files with the right file naming convention
tsgcurr = sprintf(params.TSG_FILE_FORMAT, datestr(datenum(params.YEAR, params.MONTH, params.DAY), 'yymmdd'));
tsgnext = sprintf(params.TSG_FILE_FORMAT, datestr(datenum(params.YEAR, params.MONTH, params.DAY+1), 'yymmdd'));
tsgprev = sprintf(params.TSG_FILE_FORMAT, datestr(datenum(params.YEAR, params.MONTH, params.DAY-1), 'yymmdd'));
tsgFiles = getFilesNextPrev(params.TSG_DIRECTORY, tsgcurr, tsgnext, tsgprev);

%%
% Create the TSGFileLoader object, pass in the right import method name
tfl = NAAMES_TSGFileLoader( tsgFiles, params.TSG_IMPORT_METHOD_NAME, 'temperature', ...
    params.TEMP_UNITS, 'salinity', params.SAL_UNITS, 'gps', ...
    params.GPS_UNITS );
tsgData = tfl.loadData();

%%
% loadData returns cell array of AncillaryDataobjects:
dataTypes = tfl.loadData();

%%
% loop through cell array and display some properties of each data type contained
disp('Starting looping through data returned by loadData()')

nDataTypes = length(dataTypes);

for iDataTypes = 1:nDataTypes;
   
    tempDataType = dataTypes{iDataTypes};
    fprintf('Data Type Name: %s\n', tempDataType.Name)
    fprintf('Data Type: %s\n', tempDataType.Type)
    
    tempDataType.qaqc()
    
    figure(iDataTypes)
    tempDataType.plotData()
    disp('max:')
    datestr(max(tempDataType.DataObject.Time), 'dd-mm-yyyy HH:MM:SS')
    disp('min:')
    datestr(min(tempDataType.DataObject.Time), 'dd-mm-yyyy HH:MM:SS')
    hold on
end
