%RUN_CODE - Script to run the ACS code)
%Loads configuration parameters from accode.ini and then calls the 
%following scripts:  IngestManager, PreProcessingManager, ProcessingManager
%and OutputManager
%
%
% Other m-files required: accode.ini, IngestManager, PreProcessingManager
%                         ProcessingManager, OutputManager
% Subfunctions: ini2struct
% MAT-files required: none
%
% See also: ini2struct

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% Dec 2015; Last revision: 10-12-15

%------------- BEGIN CODE --------------


% close all preexisting figures
% clear all variables and classes from workspace
close all;
% clear all;
% clear classes;

% load the configuration variables set in the .ini file
params = ini2struct('accode.ini');

% for each yearday listed in accode.ini
% process data
for iYEAR_DAY = 1:length(params.INGEST.YEAR_DAYS)
    
    params.INGEST.YEAR_DAY = params.INGEST.YEAR_DAYS(iYEAR_DAY);
    sprintf('RUN_CODE: processing year day: %u', params.INGEST.YEAR_DAY);

    %%
    % -------------------------------------------------------------------------
    % Calculate some additional DATE variables
    params.INGEST.YEAR_DAY_AS_DATENUM = doy2date( params.INGEST.YEAR_DAY, ...
        params.INGEST.YEAR);
    params.INGEST.YEAR_DAY_AS_DATEVEC = datevec( params.INGEST.YEAR_DAY_AS_DATENUM );
    params.INGEST.MONTH = params.INGEST.YEAR_DAY_AS_DATEVEC(2);
    params.INGEST.DAY = params.INGEST.YEAR_DAY_AS_DATEVEC(3);
    params.INGEST.MONTH_TEXT = datestr(params.INGEST.YEAR_DAY_AS_DATENUM, 'mmm');

    %  The file name of the log
    params.INGEST.LOG_NAME = strcat('log', num2str(params.INGEST.YEAR), '_', ...
        num2str(params.INGEST.YEAR_DAY), '.txt');
    params.INGEST.LOG_FILE = fullfile(params.RUN.LOG_DIRECTORY, params.INGEST.LOG_NAME);

    % INGEST
    if params.INGEST.AC_FILE_LOADER_TYPE
        params.INGEST.DATA_OUTPUT_FILE = strcat('acsRAW_BIN', ...
            num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY));
    else
        params.INGEST.DATA_OUTPUT_FILE = strcat('acsRAW_DAT', ...
            num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY));
    end;

    % A directory and filename for any saved matlab files from this processing
    % Saved ac data only (no ancillary data)
    params.INGEST.AC_ONLY_OUTPUT_FILE = strcat('acONLY', num2str(params.INGEST.YEAR), ...
        '_', num2str(params.INGEST.YEAR_DAY));


    % A saved 'raw' data file of ac data, plus ancillary data, for the given
    % yearday
    params.INGEST.DATA_OUTPUT_DIRECTORY = strcat(params.INGEST.DATA_OUTPUT_MAIN_DIRECTORY, ...
        num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY));

    %%
    % create log file
    L = log4m.getLogger(params.INGEST.LOG_FILE);

    % set amount of information output to screen
    if strcmp(params.RUN.LOG_SCREEN_OUTPUT_LEVEL, 'DEBUG')
        L.setCommandWindowLevel(L.DEBUG);
    elseif strcmp(params.RUN.LOG_SCREEN_OUTPUT_LEVEL, 'INFO')
        L.setCommandWindowLevel(L.INFO);
    else strcmp(params.RUN.LOG_SCREEN_OUTPUT_LEVEL, 'ERROR')
        L.setCommandWindowLevel(L.ERROR);        
    end;
 
    if strcmp(params.RUN.LOG_FILE_OUTPUT_LEVEL, 'DEBUG')
        L.setLogLevel(L.DEBUG);
    elseif strcmp(params.RUN.LOG_FILE_OUTPUT_LEVEL, 'INFO')
        L.setLogLevel(L.INFO);
    else strcmp(params.RUN.LOG_FILE_OUTPUT_LEVEL, 'ERROR')
        L.setLogLevel(L.ERROR);        
    end;
      
    if params.RUN.INGEST
        disp('------------- Running IngestManager -------------');
        run('IngestManager.m');
    else
        disp('not running ingest');
    end;
    if params.RUN.PREPROCESSING
        disp('------------- Running PreProcessingManager -------------');
        run('PreProcessingManager.m');
    else
        disp('not running ingest');
    end;
    if params.RUN.PROCESSING
        disp('------------- Running ProcessingManager -------------');
        run('ProcessingManager.m');
    end;
    if params.RUN.OUTPUT
        disp('------------- Running OutputManager -------------');
        run('OutputManager.m');
    end;

end;

%------------- END OF CODE --------------