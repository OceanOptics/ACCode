%% testACData
%
% This script tests the object ACData
% 
% 7/26/2015
% Adding functionality to test NEW findFSWTransitions.
% --------------------------------------------------------------------------
%%
% get rid of all old objects.  use delete to get rid of handle class
close all;
clear all;

params = readIngestParameters();
devFile = ACDeviceFile( params.DEVICE_FILE_LOCATION );

%%
% generate list of AC files

accurr = sprintf(params.AC_FILE_FORMAT, params.AC_SERIAL_NUMBER, params.YEAR, params.MONTH, params.DAY);
acnext = sprintf(params.AC_FILE_FORMAT, params.AC_SERIAL_NUMBER, params.YEAR, params.MONTH, params.DAY+1);
acprev = sprintf(params.AC_FILE_FORMAT, params.AC_SERIAL_NUMBER, params.YEAR, params.MONTH, params.DAY-1);
acFiles = getFilesNextPrev(params.AC_DIRECTORY, accurr, acnext, acprev);
prepacsOut = fullfile(params.PREPACS_OUTPUT_DIRECTORY, params.PREPACS_OUTPUT_FILE);

%%
% load data from AC files
%% ACFileLoader(fileNameListIn, importmethodIn, deviceFileIn, outputLocationIn)
acFileName = fullfile(params.DATA_OUTPUT_DIRECTORY, params.AC_ONLY_OUTPUT_FILE);
if params.USE_ACFILELOADER == true
    if strcmp( params.AC_FILE_LOADER_TYPE, 'bin')
        acfl = ACFileLoader( acFiles, params.AC_IMPORT_METHOD_NAME, ...
            params.DEVICE_FILE_LOCATION, prepacsOut, params.AC_UNITS, ...
            devFile.NumberWavelengths);
    else
        acfl = ACFileLoaderDat( acFiles, params.AC_IMPORT_METHOD_NAME, ...
            params.DEVICE_FILE_LOCATION, prepacsOut, params.AC_UNITS, ...
            devFile.NumberWavelengths);
    end
        acData = acfl.loadData();
    
    % if we want to save this AC data to disk, to keep from having to
    % reprocess
    if params.SAVE_AC_TO_DISK == true
        save( acFileName, 'acData');
    end
else   % not calling ACFileLoader
    acFileName = fullfile(params.DATA_OUTPUT_DIRECTORY, params.AC_ONLY_OUTPUT_FILE);
    acData = load(acFileName);
end

%%
aRawData = acData{1};
cRawData = acData{2};

% display variables of new object
aRawData.Name
aRawData.Type
aRawData.DataObject

% test methods are being called properly
figure( 1);
aRawData.plotData()
figure( 2);
cRawData.plotData();

