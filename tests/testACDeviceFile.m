% TestACDeviceFile
%
% This script tests the ACDeviceFile object

clear all;
close all;
params = readIngestParameters();

% set filename
fileNameIn = params.DEVICE_FILE_LOCATION;

            % get copy of logger
            L = log4m.getLogger();
             
          

                % check there is a file list              
                if ~isempty( fileNameIn ) 
                    FileName = fileNameIn;
                    disp(fileNameIn)
                else
                    L.error('ACDeviceFile.ACDeviceFile', 'Need file name');
                end   
                
                
                fid = fopen(fileNameIn);
                
                if fid == -1
                    L.error('ACDeviceFile.ACDeviceFile', 'Device file not successfully opened')
                end;
                      
                    %% Read first 9 lines of file
                    Intro = textscan( fid, '%s', 9, 'Delimiter', '\n');

                    % set device name
                    DeviceName = strtrim(Intro{1}{1});

                    % set serial number
                    serialNumber = Intro{1}{2};
                    SerialNumber = strtrim(strtok(serialNumber, ';'));

                    % get version number
                    versionNumber = Intro{1}{3};
                    StructureVersionNumber = strtrim(strtok(versionNumber, ';'));

                    % get calibration temps
                    cals = Intro{1}{4};
                    tempcals = textscan(cals, '"Tcal:%s C, ical: %s C. The offsets were saved to this file on %s."');
                    TCal = tempcals{1}{1};
                    ICal = tempcals{2}{1};
                    CalDate = strtok(tempcals{3}{1}, '.');

                    depth = Intro{1}{5};
                    CalibrationDepth = strtok(depth, ';');

                    baudRate = Intro{1}{6};
                    BaudRate = strtrim(strtok(baudRate, ';'));

                    pathLength = Intro{1}{7};
                    PathLength = strtrim(strtok(pathLength, ';'));

                    numWL = Intro{1}{8};
                    NumberWavelengths = strtrim(strtok(numWL, ';'));

                    numTempBins = Intro{1}{9};
                    NumberTemperatureBins = strtrim(strtok(numTempBins, ';'));

                    %%  Read Temperature Bins
                    TempBins = textscan( fid, strcat( repmat('%f', 1, str2num(NumberTemperatureBins)), ' ; temperature bins'), 1, 'Delimiter', '\n\');
                    TemperatureBins = cell2mat(TempBins(:));

                    %% Read wavelengths matrix
                    %  Size of matrix varies but follows this structure:
                    %  Field 1: Label for identifying c wavelength
                    %  Field 2: Label for identifying a wavelength
                    %  Field 3: Color for plotting within WETView (not needed)
                    %  Field 4: Clean water calibration constant for c
                    %  Field 5: Clean water calibration constant for a
                    %  The following fields contain the temperature compensation
                    %  values that correspond to the temperature bins in Line 10
                    %  The first n values are c, the next n values are for a, where
                    %  n is the variable numTempBins
                    %  There will be one line for each of the number of wavelengths
                    %  (numWL)
                    
                    numBins = str2num(NumberTemperatureBins);
                    startNextBin = numBins+1;
                    numTotalBins = (numBins)*2;
                    numRows = str2num(NumberWavelengths);
                    
                    firstcolsSpec = '%s %s %s %10.6f %10.6f';
                    endrowtext = ' "; C and A offset, and C and A temperature correction info"';
                    FormatString = strcat(firstcolsSpec, repmat('%f ', 1, numTotalBins), endrowtext);
                    
                    %% THINK ABOUT CHANGING THIS TO FSCANF IT WILL GO RIGHT INTO MATRIX
                    tempCompValueMatrix = textscan( fid, FormatString, numRows );

                    cWavelengths = cell2mat(tempCompValueMatrix{:,1});
                    aWavelengths = cell2mat(tempCompValueMatrix{:,2});
                    cWavelengths = str2num(cWavelengths(:,2:end));
                    aWavelengths = str2num(aWavelengths(:,2:end));
                    cCleanWaterCalConstant = tempCompValueMatrix{:,4};
                    aCleanWaterCalConstant = tempCompValueMatrix{:,5};
                    tempCompVals = tempCompValueMatrix(:,6:end);
                    
                    
                    
                    cTempCompensationVals = cell2mat(tempCompVals(:, 1:numBins));
                    aTempCompensationVals = cell2mat(tempCompVals(:, startNextBin:end));

                    %% Read Last line of file -- doesn't seem to contain any data
                    FooterFormatString = strcat( repmat('%f ', 1, 11), ' ; ', repmat(' %s ', 1, 11));
                    Footer = textscan( fid, FooterFormatString, 1);
                    aMaxNoise = Footer{:,1};
                    cMaxNoise = Footer{:,2}; 
                    aMaxNonConform = Footer{:,3};
                    cMaxNonConform = Footer{:,4};
                    aMaxDifference = Footer{:,5};
                    cMaxDifference = Footer{:,6};
                    aMinCounts = Footer{:,7};
                    cMinCounts = Footer{:,8};
                    rMinCounts = Footer{:,9};
                    maxTempSdev = Footer{:,10};
                    maxDepthSdev = Footer{:,11};
                    
                    %% Close file
                    closeresult = fclose(fid);
                    if closeresult == 0 
                        L.info('ACDeviceFile.ACDeviceFile', 'File close successful')
                    else
                        L.error('ACDeviceFile.ACDeviceFile', 'File close not successful')
                    end
                
%                 end   % end if file not opened



% create ACDeviceFile object
devFile = ACDeviceFile( fileNameIn );

% display info
devFile.getInfo()
