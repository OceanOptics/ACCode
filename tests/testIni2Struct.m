% close all preexisting figures and clear all variables and classes from 
% workspace
close all;
clear all;
clear classes;

% load the configuration variables set in the .ini file
params = ini2struct('accode.ini');


for iYEAR_DAY = 1:length(params.INGEST.YEAR_DAYS)
    
    params.INGEST.YEAR_DAY = params.INGEST.YEAR_DAYS(iYEAR_DAY);
    disp(params.INGEST.YEAR_DAY)


    %%
    % -------------------------------------------------------------------------
    % Calculate some additional DATE variables - DON'T EDIT
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
    % (set to L.INFO to see less)
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
        disp('going to run ingest');
        run('IngestManager.m');
    else
        disp('not running ingest');
    end;
    if params.RUN.PREPROCESSING
        disp('going to run preprocessing');
        run('PreProcessingManager.m');
    else
        disp('not running ingest');
    end;
    if params.RUN.PROCESSING
        disp('going to run processing');
        run('ProcessingManager.m');
    end;
    if params.RUN.OUTPUT
        disp('going to run output');
        run('OutputManager.m');
    end;

end;