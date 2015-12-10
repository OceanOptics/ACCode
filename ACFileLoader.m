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

%------------- BEGIN CODE --------------
classdef ACFileLoader
    
    properties
        
        FileNameList
        ImportMethodName
        DeviceFileName
        OutputLocation
        Units
        numWL
        
    end
    properties (SetAccess = private, GetAccess = private)
        L;
    end
    methods
        

        %%# ACFileLoader
        function obj = ACFileLoader(fileNameListIn, importmethodIn, ...
                deviceFileIn, outputLocationIn, unitsIn, numWLIn)
            
            %#ACFileLoader creates the data objects for the imported data
            %#
            %# SYNOPSIS ACFileLoaderDat(fileNameListIn, importmethodIn, deviceFileIn, outputLocationIn, unitsIn)
            %# INPUT fileNameListIn - a list of filenames to load
            %#    importMethodIn    - a string specifying the name of the import function
            %#    to use to read the fileNameIn
            %#    deviceFileIn      - the device file object for this data
            %#    outputLocationIn  - 
            %#    unitsIn           -
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

            for iFiles = 1:nFiles;
                
                fileName = obj.FileNameList{iFiles};
                
                obj.L.info('ACFileLoader.loadData()', ...
                    sprintf('Iteration %u.  Filename: %s', iFiles, fileName'));
               
                % -------------------------------------------------------
                % import .bin file with prepacs.exe
                
                dos(['PREPACS.EXE ' obj.DeviceFileName ' ' fileName ' ' obj.OutputLocation]);
                
                % open the data file and retrieve data from it
                fh = str2func(obj.ImportMethodName);
                origDataMatrix = fh(obj.OutputLocation);
                
                % -------------------------------------------------------            
                % Data Processing Section:
                % Logic is here because this is the object that knows about
                % the specific FILE types -- 
%%
                % get rid of any lines with a negative timestamp
                index = origDataMatrix(:,1) >= 0;
                dataMatrix = origDataMatrix(index, :);
%%                
                obj.L.debug('ACFileLoader.loadData()', ...
                    sprintf('# lines with ts less than 0: %u', size(origDataMatrix) - size(dataMatrix)));
                
                timestampsMS = dataMatrix(:,1);
                
                obj.L.debug('ACFileLoader.loadData()', ...
                    sprintf('timestampsMS start %s', datestr(timestampsMS(1), 'HH:MM:SS:FFF')));
                
                obj.L.debug('ACFileLoader.loadData()', ...
                    sprintf('timestampsMS end %s', datestr(timestampsMS(end), 'HH:MM:SS:FFF')));
         
                %% why is this hardcoded??
%                 disp('dataMatrix size:')
%                 size(dataMatrix)
                endCat = obj.numWL+1;
                cRawDataMatrix = dataMatrix(:,2:endCat);
%                 size(cRawDataMatrix)
                startAat = endCat+1;
                endAat = startAat + (obj.numWL - 1);
                aRawDataMatrix = dataMatrix(:,startAat:endAat);
%                 dataMatrix(:,endAat+1:end)
%                 size(aRawDataMatrix)
                
                % create the correct timestamp array, using the dates from the file name
                % (date not stored in timestamp)
                % and adding the timestamps from the acdata file
                [~, name, ext] = fileparts(fileName);
                shortfilename = strcat(name, ext);
                variablesFromFileName = textscan(shortfilename, 'acs%03s_%14s.bin');
                dateFromFileName = variablesFromFileName{1,2};
                fileDate = datenum(dateFromFileName, 'yyyymmddHHMMSS');
                
                timestamps = fileDate + datenum(0,0,0,0,0, (timestampsMS - timestampsMS(1))/1000);

                temp(iFiles).timestamps = timestamps;
                temp(iFiles).fileName = fileName;
                temp(iFiles).cRawDataMatrix = cRawDataMatrix;
                temp(iFiles).aRawDataMatrix = aRawDataMatrix;
            end;  % end going through list of files and loading

            % ------------------------------------------------------------
            nFiles = length(temp);
            
            obj.L.debug('ACFileLoader.loadData()', ...
                sprintf('checking timestamps.  Number of files to iterate through: %u', nFiles));

            for iFiles = 1:nFiles
                
                obj.L.info('ACFileLoader.loadData()', ...
                    sprintf('iteration %u: filename: %s', iFiles, temp(iFiles).fileName));
                
                obj.L.info('ACFileLoader.loadData()', ...
                    sprintf('----------------------------------------------------------'));

                
                [~, name, ext] = fileparts(temp(iFiles).fileName);
                shortfilename = strcat(name, ext);
                
                
                % if not last file
                if iFiles ~= nFiles
                    %if the difference between the first timestamp of the
                    %next file and the last timestamp of this file is more
                    %than 2 hours, don't append (assume next day file?)
                    
                    obj.L.debug('ACFileLoader.loadData()', 'not last file');
                    
                    if ( temp(iFiles +1).timestamps(1) - temp(iFiles).timestamps(end) ) > datenum(0,0,0,2,0,0);
                       obj.L.error('ACFileLoader.loadData()', 'greater than 2 hour gap between this and next file;');
                       continue;
                    end
                    if iFiles > 1
                        
                        % if not first file, check this file ok
                        if (temp(iFiles - 1).timestamps(end) >= temp(iFiles).timestamps(1))
                            
                            obj.L.debug('ACFileLoader.loadData()', ...
                                '***************************************************');
                            obj.L.debug('ACFileLoader.loadData()', ...
                                'this file''s first TS is BEFORE the end of the previous file''s last TS');
                            obj.L.debug('ACFileLoader.loadData()', ...
                                sprintf('this files first ts: %s', temp(iFiles).timestamps(1)));
                            obj.L.debug('ACFileLoader.loadData()', ...                            
                                sprintf('prev files last ts: %s', temp(iFiles - 1).timestamps(end)));
                            
                            if (temp(iFiles - 1).timestamps(end) > temp(iFiles).timestamps(end))
                                obj.L.error('ACFileLoader.loadData()', ...
                                    'this file''s last TS is also BEFORE the end of the previous file''s last TS');
                                continue;
                            else
                                obj.L.debug('ACFileLoader.loadData()', ...
                                    'this file''s last TS AFTER end of previous file''s last TS but there must be some overlap on beginning');

                                % find nearest last prev ts in this file, skip
                                % overlapping records
                                prevEndTime = temp(iFiles - 1).timestamps(end);
                                
                                obj.L.debug('ACFileLoader.loadData()', ...
                                    sprintf('previous end time: %s', datestr(prevEndTime)));
                                
                                [~,ind1] = min(abs(datenum(temp(iFiles).timestamps)-datenum(prevEndTime)));
                                closestStartTime = temp(iFiles).timestamps(ind1,:);
                                
                                obj.L.debug('ACFileLoader.loadData()', ...
                                    sprintf('closestStartTime: %s', datestr(closestStartTime)));
                                obj.L.debug('ACFileLoader.loadData()', ...
                                    sprintf('ind1: %u', ind1));
                                obj.L.debug('ACFileLoader.loadData()', ...
                                    sprintf('check index.  difference: %s', datestr(closestStartTime - temp(iFiles).timestamps(1))));
                                
                                % just blanking
                                temp(iFiles).timestamps = temp(iFiles).timestamps(ind1+1:end, :);
                                temp(iFiles).cRawDataMatrix = temp(iFiles).cRawDataMatrix(ind1+1:end, :);
                                temp(iFiles).aRawDataMatrix = temp(iFiles).aRawDataMatrix(ind1+1:end, :);

                            end;
                        else
                            obj.L.debug('ACFileLoader.loadData()', ...
                                sprintf('this files last timestamp should be after prev files first ts'));
                            obj.L.debug('ACFileLoader.loadData()', ...
                                sprintf('time on last: %s', datestr(temp(iFiles - 1).timestamps(end), 'HH:MM:SS:FFF')));
                            obj.L.debug('ACFileLoader.loadData()', ...
                                sprintf('time on first: %s', datestr(temp(iFiles).timestamps(1), 'HH:MM:SS:FFF')));

                        end;   % end check if last timestamp of previous file after first time stamp of this file
                    end;   % end check if not first file
                    
                    obj.L.debug('ACFileLoader.loadData()', 'outside if of check date is greater than 2 hours');
               
                end;  % if not last file

                obj.L.debug('ACFileLoader.loadData()', 'should be ok to add file');


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