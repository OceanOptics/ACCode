%ACFileLoaderDat - Opens an AC data file of type "dat" (as opposed to .bin
%file).  Iterates through a list of files specified in the 
%constructor and uses the specified import method to create new ACData
%objects.  This is the file where the logic of specifying how
%the output of the import method gets stored in the ACData object.
%
% Syntax:  loader = ACFileLoaderDat(fileNameIn, importMethodIn, deviceFileIn, outputLocationIn, unitsIn);
%
% Inputs:
%    fileNameListIn - a list of filenames to load
%    importMethodIn - a string specifying the name of the import function
%    to use to read the fileNameIn
%    deviceFileIn - the device file object for this data
%    outputLocationIn - 
%    unitsIn - 
%
% Example: 
%    acfl = ACFileLoaderDat( acFiles, params.AC_IMPORT_METHOD_NAME, params.DEVICE_FILE_LOCATION, prepacsOut, params.AC_UNITS);
%
% Other m-files required: ACData
% Subfunctions: none
% MAT-files required: none
%
% See also: ACDeviceFile,  IngestManager, ACData

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% May 2015; Last revision: 24-Nov-15
% Fixed problem with skipping data files with totally overlapped data:
% If there was a data file whose timestamps completely overlapped the 
% previous data file's timestamps, that file was being skipped.
% The next file's checking of timestamps was happening against the skipped
% file, however, and should have been happening against the last added
% file, not the last iterated file.  I added a second iterator to keep
% track of which files had been added as good,vs. which files had been
% interated through.  -WN
% June 2016 - Changed to fix discover that Compass Software was not naming 
%           - as stated in documentation (i.e. filename reflects END time
%           - of software, not START time.

%------------- BEGIN CODE --------------
%%
classdef ACFileLoaderDat
    
    properties
        FileNameList
        ImportMethodName
        DeviceFileName
        OutputLocation
        Units
        numWL
    end
    properties (SetAccess = private, GetAccess = private)
        L;   % Logger
        
    end
    methods
        
        %%# ACFileLoaderDat
        function obj = ACFileLoaderDat(fileNameListIn, importmethodIn, ...
                deviceFileIn, outputLocationIn, unitsIn, numWLIn)
            %#ACFileLoaderDat creates the data objects for the imported data
            %#
            %# SYNOPSIS ACFileLoaderDat(fileNameListIn, importmethodIn, deviceFileIn, outputLocationIn, unitsIn)
            %# INPUT fileNameListIn - a list of filenames to load
            %#    importMethodIn    - a string specifying the name of the import function
            %#    to use to read the fileNameIn
            %#    deviceFileIn      - the device file object for this data
            %#    outputLocationIn  - 
            %#    unitsIn           -
            %# OUTPUT obj
            %#
      
            
            % get copy of logger
            obj.L = log4m.getLogger();
            
            % check input arguments
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
                
            else
                error('Supply an input argument')
            end   % end if nargin > 0
            
        end   % end constructor
        
        %%# loadData
        function dataOut = loadData(obj) 
            %loadData returns an array of AncillaryData objects.  It gets
            %the list of files set when the ACFileLoader was created.  It
            %reads the .dat files with the specified import method.
            %
            % SYNOPSIS dataOut = loadData(obj) 
            % INPUT  obj     - this object
            % OUTPUT dataOut - cell array of adata and cdata data objects
          
            nFiles = length(obj.FileNameList);
            obj.L.info('ACFileLoaderDat', ...
                sprintf('Number of files to iterate through: %u', nFiles));
            
            % iterate through files, load from original data files and put
            % in temporary structure to hold for looping & checking
            for iFiles = 1:nFiles;
                
                fileName = obj.FileNameList{iFiles};
                obj.L.info('ACFileLoaderDat', ...
                    sprintf('Iteration %u.  Filename: %s', iFiles, fileName'));
                % -------------------------------------------------------
               
                % open the data file and retrieve data from it
                fh = str2func(obj.ImportMethodName);
                disp(obj.ImportMethodName)
                %changed 11/23
                % get rid of any lines with a negative timestamp
                origDataMatrix = fh(fileName);                
                index = origDataMatrix(:,1) >= 0;
                dataMatrix = origDataMatrix(index, :);
                % end change
                
                % -------------------------------------------------------            
                % Data Processing Section:
                % Logic is here because this is the object that knows about
                % the specific FILE types -- 
                
                timestampsMS = dataMatrix(:,1);
                
                obj.L.debug('ACFileLoaderDat', ...
                    sprintf('dataMatrix size: %u x %u', size(dataMatrix)));
                
                % calculate dimensions of matrix
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
                variablesFromFileName = textscan(shortfilename, 'acs%03s_%14s.dat');
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
                timestamps = fileDate - datenum(0,0,0,0,0, (timestampsMS(end) - timestampsMS)/1000);
                %------------------------------------------------------

                temp(iFiles).timestamps = timestamps;
                temp(iFiles).fileName = fileName;
                temp(iFiles).cRawDataMatrix = cRawDataMatrix;
                temp(iFiles).aRawDataMatrix = aRawDataMatrix;
                
                obj.L.debug('ACFileLoaderDat', ...
                    sprintf('timestamps start %s', datestr(timestamps(1), 'HH:MM:SS:FFF')));
                obj.L.debug('ACFileLoaderDat', ...
                    sprintf('timestampsMS end %s', datestr(timestamps(end), 'HH:MM:SS:FFF')));

            end; %# for iFiles = 1:nFiles;

            nFiles = length(temp);
            
            obj.L.info('ACFileLoaderDat', ...
                sprintf('checking timestamps.  Number of files to iterate through: %u', nFiles));
            
            % create an index for the last good file.  skipped files will
            % not iterate this index
            iLastAddedFile = 0;            
            
            for iFiles = 1:nFiles

                currFile = temp(iFiles);

                obj.L.debug('ACFileLoaderDat', ...
                    sprintf('iteration %u: filename: %s', iFiles, currFile.fileName));
                if iLastAddedFile > 0
                    obj.L.debug('ACFileLoaderDat', ...
                        sprintf('last good added filename: %s', temp(iLastAddedFile).fileName));
                end
                [~, name, ext] = fileparts(currFile.fileName);
                shortfilename = strcat(name, ext);
                
                
                % if not last file
                if iFiles ~= nFiles
                    nextFile = temp(iFiles +1);
                    
                    obj.L.debug('ACFileLoaderDat', 'not last file');
                    
                    %if the difference between the first timestamp of the
                    %next file and the last timestamp of this file is more
                    %than 2 hours, don't append (assume next day file?)

                    if ( nextFile.timestamps(1) - currFile.timestamps(end) ) > datenum(0,0,0,2,0,0);
                        obj.L.debug('ACFileLoaderDat', ...
                            'greater than 2 hour gap between this and next file; SHOULD SKIP TO NEXT FILE');
                        continue;
                        obj.L.error('ACFileLoaderDat', ...
                            'DID NOT SKIP TO NEXT FILE IN LOOP');
                    end
                else
                    obj.L.debug('ACFileLoaderDat', 'last file');
                end  %# if not last file;
                
                if iFiles > 1
                    % if not first file, check this file ok
                    prevAddedFile = temp(iLastAddedFile);

                    %if (prevFile.timestamps(end) > currFile.timestamps(1))
                    if (currFile.timestamps(1) < prevAddedFile.timestamps(end))    

                        obj.L.debug('ACFileLoaderDat', ...
                            'this file''s first TS is BEFORE the end of the previous file''s last TS');

                        obj.L.debug('ACFileLoaderDat', ...
                            sprintf('this files first ts: %s', ...
                            datestr((currFile.timestamps(1)), 'HH:MM:SS:FFF')));
                        obj.L.debug('ACFileLoaderDat', ...
                            sprintf('prev files last ts: %s', ...
                            datestr((prevAddedFile.timestamps(end)), 'HH:MM:SS:FFF')));

                        %if (prevFile.timestamps(end) > currFile.timestamps(end))
                        if (currFile.timestamps(end) < prevAddedFile.timestamps(end))
                            obj.L.error('ACFileLoaderDat', ...
                                'curFile last TS BEFORE prevFile last TS - SKIPPING');
                            continue;
                            obj.L.error('ACFileLoaderDat', ...
                                'DID NOT SKIP TO NEXT FILE IN LOOP');
                        else

                            obj.L.error('ACFileLoaderDat', ...
                                'currFile last TS AFTER the prevFile last TS - BUT OVERLAP AT BEGINNING')

                            % find nearest last prev ts in this file, skip
                            % overlapping records

                            prevEndTime = prevAddedFile.timestamps(end);

                            obj.L.debug('ACFileLoaderDat', ...
                                sprintf('prevEndTime: %s:', datestr(prevEndTime, 'HH:MM:SS:FFF')));

                            [~,ind1] = min(abs(datenum(currFile.timestamps)-datenum(prevEndTime)));
                            closestStartTime = currFile.timestamps(ind1,:);

                            obj.L.debug('ACFileLoaderDat', ...
                                sprintf('closest start time: %s', datestr(closestStartTime, 'HH:MM:SS:FFF')));


                            obj.L.debug('ACFileLoaderDat', ...
                                sprintf('length of timestamps: %u', length(currFile.timestamps)));

                            % let's try just blanking
                            currFile.timestamps = currFile.timestamps(ind1+1:end, :);

                            obj.L.debug('ACFileLoaderDat', ...
                                sprintf('adjusted timestamps length: %u', length(currFile.timestamps)));

                            % BLANK C DATA
                            currFile.cRawDataMatrix = currFile.cRawDataMatrix(ind1+1:end, :);

                            obj.L.debug('ACFileLoaderDat', ...
                                sprintf('original A data size: %u x %u', size(currFile.aRawDataMatrix)));

                            % BLANK A DATA
                            currFile.aRawDataMatrix = currFile.aRawDataMatrix(ind1+1:end, :);

                            obj.L.debug('ACFileLoaderDat', ...
                                sprintf('adjusted A data size %u x %u', size(currFile.aRawDataMatrix)));

                            obj.L.debug('ACFileLoaderDat', 'finished adjusting for time');

                        end; %#(prevFile.timestamps(end) > currFile.timestamps(end))
                    else
                        obj.L.debug('ACFileLoaderDat', 'currFile 1st ts NOT before prevFile last ts');
                    end; %#(prevFile.timestamps(end) > currFile.timestamps(1))
                else
                    obj.L.debug('ACFileLoaderDat', 'First file');
                end;   % end check if not first file

                obj.L.debug('ACFileLoaderDat', 'outside of check date is greater than 2 hours');
                
                
                obj.L.debug('ACFileLoaderDat','should be ok to add');
                % create new ACData, constructor: ACData(nameIn, dataValuesIn, timestampsIn
                % acdata = ACData('testAC', [cRawData, aRawData], timestamps);                
                % ------------------------------------------------------
                % if a or c Data already exists
                if ismember({'aRawData'}, who)
                    obj.L.info('ACFileLoaderDat', sprintf('adding %s', shortfilename));
                    
                    %append data
                    aRawData.addData(['a: ' shortfilename], currFile.aRawDataMatrix, currFile.timestamps);
                    cRawData.addData(['c: ' shortfilename], currFile.cRawDataMatrix, currFile.timestamps);
                    disp('****')
                    iLastAddedFile = iFiles;

                else   % a or c Data doesn't already exist
                    obj.L.info('ACFileLoaderDat', sprintf('creating AC data with %s', shortfilename));              
                    aRawData = ACData(['a: ' shortfilename], currFile.aRawDataMatrix, currFile.timestamps, obj.Units);
                    cRawData = ACData(['c: ' shortfilename], currFile.cRawDataMatrix, currFile.timestamps, obj.Units);
                    iLastAddedFile = iFiles;
                    
                end   % end check for existing data object


            end % end for loop

            % create a cell array to hold both
            dataOut = {aRawData, cRawData};
            
        end   % end loadData function
    end   % end methods block
end   % end classDef
