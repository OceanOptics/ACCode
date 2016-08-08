%ACFileLoader - Open an AC data file and create two new ACData objects.
% Iterates through a list of files specified in the constructor and uses the
% specified import method to create new ACData objects.
% This is the file where the logic of specifying how
% the output of the import method gets stored in the ACData object.
%
% Syntax: loader = ACFileLoader(fileNameListIn, importmethodIn,
% deviceFileIn, outputLocationIn, unitsIn);
%
% Inputs:
%    fileNameListIn - a list of full file names of the ac files
%    importMethodIn - the name of the import script to use to read files
%    deviceFileIn - the device file object for this data
%    outputLocationIn - 
%    unitsIn - 
%
% Example:
%     acfl = ACFileLoader( acFiles, params.AC_IMPORT_METHOD_NAME, params.DEVICE_FILE_LOCATION, prepacsOut, params.AC_UNITS);
%
% Other m-files required: ACData
% Subfunctions: none
% MAT-files required: none
%
% See also: ACDeviceFile,  IngestManager, ACData
%
% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% May 2015; Last revision: 13-08-15
% May 2016 - Added ability to run from OSX
% June 2016 - Changed to fix discover that Compass Software was not naming 
%           - as stated in documentation (i.e. filename reflects END time
%           - of software, not START time.

%------------- BEGIN CODE --------------
classdef ACFileLoader
    
    properties
        
        FileNameList
        ImportMethodName
        DeviceFileName
        PrepacsBin
        OutputLocation
        Units
        numWL
        PrepacsLocation
        
    end
    properties (SetAccess = private, GetAccess = private)
        L;
    end
    methods
        

        %%# ACFileLoader
        function obj = ACFileLoader(fileNameListIn, importmethodIn, ...
                deviceFileIn, outputLocationIn, unitsIn, numWLIn, prepacsBin)
           
            %#ACFileLoader creates the data objects for the imported data
            %#
            %# SYNOPSIS ACFileLoaderDat(fileNameListIn, importmethodIn, deviceFileIn, outputLocationIn, unitsIn)
            %# INPUT fileNameListIn - a list of filenames to load
            %#    importMethodIn    - a string specifying the name of the import function
            %#    to use to read the fileNameIn
            %#    deviceFileIn      - the device file object for this data
            %#    outputLocationIn  - 
            %#    unitsIn           -
            %#    prepacsBin        - name of binary to run prepacs
            %# OUTPUT obj - this object
            %#          
            if nargin > 0
                
                % check there is a file list              
                if ~isempty( fileNameListIn ) 
                    obj.FileNameList = fileNameListIn;
                else
                    error('Need file list')
                end
                
                % check there is a import method specified
                if ischar( importmethodIn )
                    obj.ImportMethodName = importmethodIn;
                else
                    error('Need import method')
                end 
                
                % check there is a device file specified
                if ischar( deviceFileIn )
                    obj.DeviceFileName = deviceFileIn;
                else
                    error('Need device file')
                end
                
                if ischar( outputLocationIn )
                    obj.OutputLocation = outputLocationIn;
                else
                    error('Need output location for prepacs.tmp')
                end
                
                if ischar( unitsIn )
                    obj.Units = unitsIn;
                else
                    error('Need units')
                end
                obj.numWL = str2num(numWLIn);
                
                if ischar( prepacsBin )
                    obj.PrepacsBin = prepacsBin;
                else
                    error('Need prepacs binary name')
                end
            else
                error('Supply an input argument')
            end   % end if nargin > 0
            
            % create logger
            obj.L = log4m.getLogger();
            obj.L.debug('ACFileLoader.ACFileLoader()','Created object');
            
        end   % end constructor
        
        %%# loadData
        function dataOut = loadData(obj) 
            % loadData returns an array of AncillaryData objects
            % Example:  one "ac file" contains both aData and cData
            % Get the list of files set when the ACFileLoader was created.
            % Iterate through the list and call prepacs.exe to read the
            % .bin files, creating a prepacs.tmp file.
            % Read the prepacs.tmp file with the specified importfileACS.m
            % method.   
            %loadData returns an array of AncillaryData objects.  It gets
            %the list of files set when the ACFileLoader was created.  It
            %reads the .dat files with the specified import method.
            %
            % SYNOPSIS dataOut = loadData(obj) 
            % INPUT  obj     - this object
            % OUTPUT dataOut - cell array of adata and cdata data objects
            
            nFiles = length(obj.FileNameList);
            obj.L.debug('ACFileLoader.loadData()', ...
                sprintf('Number of files to iterate through: %u', nFiles));
            
            nFilesWithGoodData = 0;

            for iFiles = 1:nFiles;
                
                fileName = obj.FileNameList{iFiles};
                
                obj.L.info('ACFileLoader.loadData()', ...
                    sprintf('Iteration %u.  Filename: %s', iFiles, fileName'));
               
                % -------------------------------------------------------
                % import .bin file with prepacs.exe
                
                prepacs_output_filename = strcat(obj.OutputLocation, 'prepacs', num2str(iFiles), '.tmp');
                
                dos([obj.PrepacsBin ' ' obj.DeviceFileName ' ' fileName ' ' prepacs_output_filename]);
                
                % open the data file and retrieve data from it
                fh = str2func(obj.ImportMethodName);
                origDataMatrix = fh(prepacs_output_filename);
                
                % -------------------------------------------------------            
                % Data Processing Section:
                % Logic is here because this is the object that knows about
                % the specific FILE types -- 
                
                if ~isempty(origDataMatrix)
                    
                    % increase iterator
                    nFilesWithGoodData = nFilesWithGoodData + 1;
                    
                    % get rid of any lines with a negative timestamp
                    index = origDataMatrix(:,1) >= 0;
                    dataMatrix = origDataMatrix(index, :);
   
                    obj.L.debug('ACFileLoader.loadData()', ...
                        sprintf('# lines with ts less than 0: %u', size(origDataMatrix) - size(dataMatrix)));

                    timestampsMS = dataMatrix(:,1);

                    endCat = obj.numWL+1;
                    cRawDataMatrix = dataMatrix(:,2:endCat);
                    startAat = endCat+1;
                    endAat = startAat + (obj.numWL - 1);
                    aRawDataMatrix = dataMatrix(:,startAat:endAat);

                    % create the correct timestamp array, using the dates from the file name
                    % (date not stored in timestamp)
                    % and adding the timestamps from the acdata file
                    [~, name, ext] = fileparts(fileName);
                    shortfilename = strcat(name, ext);
                    variablesFromFileName = textscan(shortfilename, 'acs%03s_%14s.bin');
                    dateFromFileName = variablesFromFileName{1,2};
                    fileDate = datenum(dateFromFileName, 'yyyymmddHHMMSS');

                    %------------------------------------------------------
                    %22Jun16:
                    % changed code to reflect issue with compass software
                    % realized filenames are not named according to the 
                    % beginning time of the datafile name, but rather, 
                    % the END of times in data file
                    
                    % OLD CODE:
                    %timestamps = fileDate + datenum(0,0,0,0,0, (timestampsMS - timestampsMS(1))/1000);

                    % NEW CODE:
                    timestamps = fileDate - datenum(0,0,0,0,0, ...
                        (timestampsMS(end) - timestampsMS)/1000);
                    %------------------------------------------------------
                    temp(nFilesWithGoodData).timestamps = timestamps;
                    temp(nFilesWithGoodData).fileName = fileName;
                    temp(nFilesWithGoodData).cRawDataMatrix = cRawDataMatrix;
                    temp(nFilesWithGoodData).aRawDataMatrix = aRawDataMatrix;
                else
                    % don't increase iterator, not good file
                    obj.L.debug('ACFileLoader.loadData()', ...
                        sprintf('No data: %u', nFiles));
                end;
                
            end;  % for iFiles = 1:nFiles -- end going through list of files and loading

            % ------------------------------------------------------------
            nGoodFiles = length(temp);
            numFilesAdded = 0;
           
            lastFileAddedFirstTS = NaN;
            lastFileAddedLastTS = NaN;
            obj.L.debug('ACFileLoader.loadData()', ...
                sprintf('checking timestamps.  Number of files to iterate through: %u', nGoodFiles));

            for iFiles = 1:nGoodFiles
                [~, name, ext] = fileparts(temp(iFiles).fileName);
                shortfilename = strcat(name, ext);
                
                obj.L.info('ACFileLoader.loadData()', ...
                    sprintf('iteration %u: filename: %s', iFiles, shortfilename));
                obj.L.info('ACFileLoader.loadData()', ...
                    sprintf('----------------------------------------------------------'));
                % set TS vars
                thisFilesFirstTS = temp(iFiles).timestamps(1);                
                thisFilesFirstTS_DV = datevec(temp(iFiles).timestamps(1));
                
                thisFilesLastTS = temp(iFiles).timestamps(end);
                thisFilesLastTS_DV = datevec(temp(iFiles).timestamps(end));
                
                % if not first file and a file has already been added
                if (1 < iFiles) && (iFiles <= nGoodFiles) && numFilesAdded >= 1 

                    if (lastFileAddedLastTS > thisFilesFirstTS)
                         obj.L.debug('ACFileLoader.loadData()', ...
                            '***************************************************');
                        obj.L.debug('ACFileLoader.loadData()', ...
                            'this file''s first TS is BEFORE the end of the previous file''s last TS.');
                        obj.L.debug('ACFileLoader.loadData()', ...
                            sprintf('this files first ts: %s', datestr(thisFilesFirstTS)));  %datestr(temp(iFiles).timestamps(1))));
                        obj.L.debug('ACFileLoader.loadData()', ...                            
                            sprintf('prev files last ts: %s', datestr(lastFileAddedLastTS))); %temp(iFiles - 1).timestamps(end))));
                        if lastFileAddedLastTS > thisFilesLastTS
                            obj.L.error('ACFileLoader.loadData()', ...
                                'ERROR: this file''s last TS is also BEFORE the end of the previous file''s last TS');
                            % don't add
                            continue;
                        else
                            obj.L.debug('ACFileLoader.loadData()', ...
                                'this file''s last TS AFTER end of last added files''s last TS but there must be some overlap on beginning');

                            % find nearest last prev ts in this file, skip
                            % overlapping records
                            obj.L.debug('ACFileLoader.loadData()', ...
                                sprintf('previous end time: %s', datestr(lastFileAddedLastTS)));
                                
                            [~,ind1] = min(abs(datenum(temp(iFiles).timestamps)-datenum(lastFileAddedLastTS)));
                            closestStartTime = temp(iFiles).timestamps(ind1,:);
                                
                            obj.L.debug('ACFileLoader.loadData()', ...
                                sprintf('closestStartTime: %s', datestr(closestStartTime)));
                            obj.L.debug('ACFileLoader.loadData()', ...
                                sprintf('ind1: %u', ind1));
                            obj.L.debug('ACFileLoader.loadData()', ...
                                sprintf('check index.  difference: %s', datestr(closestStartTime - thisFilesFirstTS)));
                                
                            % just blanking
                            temp(iFiles).timestamps = temp(iFiles).timestamps(ind1+1:end, :);
                            temp(iFiles).cRawDataMatrix = temp(iFiles).cRawDataMatrix(ind1+1:end, :);
                            temp(iFiles).aRawDataMatrix = temp(iFiles).aRawDataMatrix(ind1+1:end, :);

                            end;  %if lastFileAddedLastTimestamp > thisFileLastTS
                    else %(lastFileAddedLastTimestamp > thisFilesFirstTS)
                        obj.L.debug('ACFileLoader.loadData()', ...
                            sprintf('this files last timestamp is AFTER prev files first ts'));
                        obj.L.debug('ACFileLoader.loadData()', ...
                            sprintf('prev file last ts: %s', ...
                            datestr(lastFileAddedLastTS, 'YY-mm-DD HH:MM:SS:FFF')));
                        obj.L.debug('ACFileLoader.loadData()', ...
                            sprintf('this file first ts: %s', ...
                            datestr(thisFilesFirstTS, 'YY-mm-DD HH:MM:SS:FFF')));
                    end;   % end check if last timestamp of previous file after first time stamp of this file
                    if iFiles == nGoodFiles  % last file to be added
                        % check this last file isn't more than 2 hours past
                        % previous
                        obj.L.debug('ACFileLoader.loadData()', 'Last file');
                        lastFileAddedLastTS_DV = datevec(lastFileAddedLastTS);
                        if ((thisFilesFirstTS - lastFileAddedLastTS) > datenum(0,0,0,2,0,0)) ...
                                && (thisFilesFirstTS_DV(3) ~= lastFileAddedLastTS_DV(3))
                           obj.L.error('ACFileLoader.loadData()', ...
                               'ERROR: 2+ hr gap between last file and last added, Also: DIFFERENT DAY;');
                        else
                            obj.L.debug('ACFileLoader.loadData()', ...
                                'Last file within 2 hours of last added file');
                        end;
                    end;


                elseif numFilesAdded == 0;   %if this is the first file to be added
                    obj.L.debug('ACFileLoader.loadData()', 'No files added yet');
                    
                    %if the difference between the first timestamp of the
                    %next file and the last timestamp of this file is more
                    %than 2 hours, and is 100% a different day, don't use
                    
                    if size(temp,2) > 1
                      nextFilesFirstTimestampDV = datevec(temp(iFiles+1).timestamps(1));
                      nextFilesFirstTimestamp = temp(iFiles+1).timestamps(1);

                      if ((nextFilesFirstTimestamp - thisFilesLastTS) > datenum(0,0,0,2,0,0)) ...
                              && nextFilesFirstTimestampDV(3) ~= thisFilesLastTS_DV(3);
                         obj.L.error('ACFileLoader.loadData()', 'ERROR: greater than 2 hour gap between this and next file, and DIFFERENT DAY;');
                         continue;
                      else
                          obj.L.debug('ACFileLoader.loadData()', 'First file within 2 hours of next OR same day');
                      end;
                    else
                      obj.L.debug('ACFileLoader.loadData()', 'Only one file in temp');
                    end;
                end;  % (1 < iFiles <= nGoodFiles)  && numFilesAdded >= 1 

                obj.L.debug('ACFileLoader.loadData()', 'should be ok to add file');
                
                % update 'lastFileAdded' timestamps
                lastFileAddedLastTS = thisFilesLastTS;
                numFilesAdded = numFilesAdded + 1;
                
                % create new ACData, constructor: ACData(nameIn, dataValuesIn, timestampsIn
                % acdata = ACData('testAC', [cRawData, aRawData], timestamps);                
                % ------------------------------------------------------
                % if a or c Data already exists
                if ismember({'aRawData'}, who)
                    obj.L.info('ACFileLoader.loadData()', ...
                        sprintf('adding %s', shortfilename));
                    
                    %append data
                    aRawData.addData(['a: ' shortfilename], temp(iFiles).aRawDataMatrix, temp(iFiles).timestamps);

                    cRawData.addData(['c: ' shortfilename], temp(iFiles).cRawDataMatrix, temp(iFiles).timestamps);

                else   % a or c Data doesn't already exist
                    obj.L.info('ACFileLoader.loadData()', ...
                         sprintf('creating AC data with %s', shortfilename) ); 

                    aRawData = ACData(['a: ' shortfilename], temp(iFiles).aRawDataMatrix, temp(iFiles).timestamps, obj.Units);
                    cRawData = ACData(['c: ' shortfilename], temp(iFiles).cRawDataMatrix, temp(iFiles).timestamps, obj.Units);

                end   % end check for existing data object

            end % end for loop

            % create a cell array to hold both
            dataOut = {aRawData, cRawData};
            obj.L.info('ACFileLoader.loadData()','End Method');
            
        end   % end loadData function
    end   % end methods block
end   % end classDef

%------------- END CODE --------------