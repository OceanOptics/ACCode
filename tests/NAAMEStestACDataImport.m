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
L = log4m.getLogger();
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

            for iFiles = 1:length(acFiles);
                
                fileName = acFiles{iFiles};
                
                L.info('ACFileLoader.loadData()', ...
                    sprintf('Iteration %u.  Filename: %s', iFiles, fileName'));
               
                % -------------------------------------------------------
                % import .bin file with prepacs.exe
                iFiles
                numstr = num2str(iFiles)
                tempfilename = strcat('prepacs',numstr,'.tmp')
                tempfullfilename = strcat('C:\Users\Wendy\Documents\data\NAAMES\TEMP\', tempfilename)
                dos(['PREPACS.EXE ' params.DEVICE_FILE_LOCATION ' ' fileName ' ' tempfullfilename]);
                
            end;
     
            %%
            for iFiles = 1:length(acFiles);
                
                % open the data file and retrieve data from it
                fh = str2func('importfileACS');
                empfilename = strcat('prepacs',num2str(i),'.tmp')
                tempfullfilename = strcat('C:\Users\Wendy\Documents\data\NAAMES\TEMP\', tempfilename)
                origDataMatrix = fh(tempfullfilename);
                
                
                
                % -------------------------------------------------------            
                % Data Processing Section:
                % Logic is here because this is the object that knows about
                % the specific FILE types -- 
%
                % get rid of any lines with a negative timestamp
                index = origDataMatrix(:,1) >= 0;
                dataMatrix = origDataMatrix(index, :);
%                
                L.debug('ACFileLoader.loadData()', ...
                    sprintf('# lines with ts less than 0: %u', size(origDataMatrix) - size(dataMatrix)));
                
                timestampsMS = dataMatrix(:,1);
                
                L.debug('ACFileLoader.loadData()', ...
                    sprintf('timestampsMS start %s', datestr(timestampsMS(1), 'HH:MM:SS:FFF')));
                
                L.debug('ACFileLoader.loadData()', ...
                    sprintf('timestampsMS end %s', datestr(timestampsMS(end), 'HH:MM:SS:FFF')));
         
                %% why is this hardcoded??
                disp('dataMatrix size:')
                size(dataMatrix)
                endCat = 77
                cRawDataMatrix = dataMatrix(:,2:endCat);
                size(cRawDataMatrix)
                startAat = endCat+1
                endAat = startAat + (77 - 1)
                aRawDataMatrix = dataMatrix(:,startAat:endAat);
%                 dataMatrix(:,endAat+1:end)
                size(aRawDataMatrix)
                
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
%%
            % start from aCFeilloader
            nFiles = length(temp);
            
            L.debug('ACFileLoader.loadData()', ...
                sprintf('checking timestamps.  Number of files to iterate through: %u', nFiles));

            for iFiles = 1:nFiles
                
                L.info('ACFileLoader.loadData()', ...
                    sprintf('iteration %u: filename: %s', iFiles, temp(iFiles).fileName));
                
                L.info('ACFileLoader.loadData()', ...
                    sprintf('----------------------------------------------------------'));

                temp(iFiles)
                [~, name, ext] = fileparts(temp(iFiles).fileName);
                shortfilename = strcat(name, ext);
                
                
                % if not last file
                if iFiles ~= nFiles
                    %if the difference between the first timestamp of the
                    %next file and the last timestamp of this file is more
                    %than 2 hours, don't append (assume next day file?)
                    
                    L.debug('ACFileLoader.loadData()', 'not last file');
                    
                    if ( temp(iFiles +1).timestamps(1) - temp(iFiles).timestamps(end) ) > datenum(0,0,0,2,0,0);
                       L.error('ACFileLoader.loadData()', 'greater than 2 hour gap between this and next file;');
                       continue;
                    end
                    if iFiles > 1
                        
                        % if not first file, check this file ok
                        if (temp(iFiles - 1).timestamps(end) >= temp(iFiles).timestamps(1))
                            
                            L.debug('ACFileLoader.loadData()', ...
                                '***************************************************');
                            L.debug('ACFileLoader.loadData()', ...
                                'this file''s first TS is BEFORE the end of the previous file''s last TS');
                            L.debug('ACFileLoader.loadData()', ...
                                sprintf('this files first ts: %s', temp(iFiles).timestamps(1)));
                            L.debug('ACFileLoader.loadData()', ...                            
                                sprintf('prev files last ts: %s', temp(iFiles - 1).timestamps(end)));
                            
                            if (temp(iFiles - 1).timestamps(end) > temp(iFiles).timestamps(end))
                                L.error('ACFileLoader.loadData()', ...
                                    'this file''s last TS is also BEFORE the end of the previous file''s last TS');
                                continue;
                            else
                                L.debug('ACFileLoader.loadData()', ...
                                    'this file''s last TS AFTER end of previous file''s last TS but there must be some overlap on beginning');

                                % find nearest last prev ts in this file, skip
                                % overlapping records
                                prevEndTime = temp(iFiles - 1).timestamps(end);
                                
                                L.debug('ACFileLoader.loadData()', ...
                                    sprintf('previous end time: %s', datestr(prevEndTime)));
                                
                                [~,ind1] = min(abs(datenum(temp(iFiles).timestamps)-datenum(prevEndTime)));
                                closestStartTime = temp(iFiles).timestamps(ind1,:);
                                
                                L.debug('ACFileLoader.loadData()', ...
                                    sprintf('closestStartTime: %s', datestr(closestStartTime)));
                                L.debug('ACFileLoader.loadData()', ...
                                    sprintf('ind1: %u', ind1));
                                L.debug('ACFileLoader.loadData()', ...
                                    sprintf('check index.  difference: %s', datestr(closestStartTime - temp(iFiles).timestamps(1))));
                                
                                % just blanking
                                temp(iFiles).timestamps = temp(iFiles).timestamps(ind1+1:end, :);
                                temp(iFiles).cRawDataMatrix = temp(iFiles).cRawDataMatrix(ind1+1:end, :);
                                temp(iFiles).aRawDataMatrix = temp(iFiles).aRawDataMatrix(ind1+1:end, :);

                            end;
                        else
                            L.debug('ACFileLoader.loadData()', ...
                                sprintf('this files last timestamp should be after prev files first ts'));
                            L.debug('ACFileLoader.loadData()', ...
                                sprintf('time on last: %s', datestr(temp(iFiles - 1).timestamps(end), 'HH:MM:SS:FFF')));
                            L.debug('ACFileLoader.loadData()', ...
                                sprintf('time on first: %s', datestr(temp(iFiles).timestamps(1), 'HH:MM:SS:FFF')));

                        end;   % end check if last timestamp of previous file after first time stamp of this file
                    end;   % end check if not first file
                    
                    L.debug('ACFileLoader.loadData()', 'outside if of check date is greater than 2 hours');
               
                end;  % if not last file

                L.debug('ACFileLoader.loadData()', 'should be ok to add file');


                % create new ACData, constructor: ACData(nameIn, dataValuesIn, timestampsIn
                % acdata = ACData('testAC', [cRawData, aRawData], timestamps);                
                % ------------------------------------------------------
                % if a or c Data already exists
                if ismember({'aRawData'}, who)
                    L.info('ACFileLoader.loadData()', ...
                        sprintf('adding %s', shortfilename));
                    
                    %append data
                    aRawData.addData(['a: ' shortfilename], temp(iFiles).aRawDataMatrix, temp(iFiles).timestamps);

                    cRawData.addData(['c: ' shortfilename], temp(iFiles).cRawDataMatrix, temp(iFiles).timestamps);

                else   % a or c Data doesn't already exist
                    L.info('ACFileLoader.loadData()', ...
                         sprintf('creating AC data with %s', shortfilename) ); 

                    aRawData = ACData(['a: ' shortfilename], temp(iFiles).aRawDataMatrix, temp(iFiles).timestamps, 'units');
                    cRawData = ACData(['c: ' shortfilename], temp(iFiles).cRawDataMatrix, temp(iFiles).timestamps, 'units');

                end   % end check for existing data object

            end % end for loop

            % create a cell array to hold both
            dataOut = {aRawData, cRawData};
            L.info('ACFileLoader.loadData()','End Method');
            
            

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

