function [fileList] = getFilesNextPrev( directoryIn, currentSpecIn, nextSpecIn, prevSpecIn)
%GETFILESNEXTPREV - Outputs a list of files available for a specified
%yearday, including the last file from the previous day and the first file
%from the next day.  These extra files are included to make sure entire
%filter periods are captured.
%
% Syntax:  [fileList] = getFilesNextPrev( directoryIn, currentSpecIn, nextSpecIn, prevSpecIn)
%
% Inputs:
%    directoryIn - Description
%    currentSpecIn - Description
%    nextSpecIn - Description
%    prevSpecIn - 
%
% Outputs:
%    fileList - 
%
% Example: 
%    flowFiles = getFilesNextPrev(params.FLOW_DIRECTORY, flowcurr, flownext, flowprev);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: INGEST_MANAGER.m

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% May 2015; Last revision: 13-08-15

%------------- BEGIN CODE --------------
    % get handle to logger
    L = log4m.getLogger();
    
    % need check for inputs here
    if nargin == 4
                
        validateattributes( directoryIn, {'char'}, {'nonempty'} );
        validateattributes( currentSpecIn, {'char'}, {'nonempty'} );
        validateattributes( nextSpecIn, {'char'}, {'nonempty'} );
        validateattributes( prevSpecIn, {'char'}, {'nonempty'} );
                
    else
        L.error('getFilesNextPrev', 'function requires 4 inputs');
    end;

    % generate list of files
    currentDayFiles = dir( fullfile( directoryIn, currentSpecIn));
    currentDayFiles = {currentDayFiles.name}';
    currentDayFiles = fullfile(directoryIn, currentDayFiles);

    for i = 1:length(currentDayFiles)
        L.info('getFilesNextPrev',  currentDayFiles{i});
    end
    
    % get file from previous day
    prevDayFiles = dir( fullfile( directoryIn, prevSpecIn));
    if isempty( prevDayFiles )
        % no files from previous day available
        L.info('getFilesNextPrev', 'no file from previous day');
    else
        % files from previous day DO exist
        % take last one
        prevDayFiles = {prevDayFiles.name}';
        prevFile = fullfile(directoryIn, prevDayFiles(end));
        L.info('getFilesNextPrev', sprintf('File from previous day exists: %s', prevFile{1}));
                
    end

    nextDayFiles = dir( fullfile( directoryIn, nextSpecIn));
    if isempty( nextDayFiles )
        % no files from next day available
        L.info('getFilesNextPrev', 'no file from next day');
    else
        % files from next day are available
        % take first one
        nextDayFiles = {nextDayFiles.name}';
        nextFile = fullfile(directoryIn,nextDayFiles(1));
        L.info('getFilesNextPrev', sprintf('File from next day exists: %s', nextFile{1}));
        
    end

    %create list of files to load
    if (ismember('prevFile', who))
        if (ismember('nextFile', who))
            % both day before and day after
            L.info('getFilesNextPrev', 'Files for both day before and day after available');
            fileList = [prevFile; currentDayFiles; nextFile];
        else
            % day before available, but not day after
            L.info('getFilesNextPrev', 'Files for day before, but not day after available');
            fileList = [prevFile; currentDayFiles];
        end
    else
        %day before not available
        if (ismember('nextFile', who))
            %day after available
            L.info('getFilesNextPrev', 'File for day after, but not day before available');
            fileList = [ currentDayFiles; nextFile];
        else
            %neither available
            L.info('getFilesNextPrev', 'Neither files for day before or day after available');
            fileList = currentDayFiles;
        end
    end
end
%------------- END OF CODE --------------