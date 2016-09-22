%ACDEVICEFILE - Object represents the data contained in an ACS Device File
%
% Syntax:  obj = ACDeviceFile(fileNameIn)
%
% Inputs:
%    fileNameIn - the full file name of the device file
%
%
% Example: 
%    devFile = ACDeviceFile( params.DEVICE_FILE_LOCATION );
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: IngestManager
%
% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% May 2015; Last revision: 13-08-15

%------------- BEGIN CODE --------------
%%
classdef ACDeviceFile < handle
    
    properties
       FileName
       FileType
       DeviceName
       SerialNumber
       StructureVersionNumber
       TCal
       ICal
       CalDate
       CalibrationDepth
       BaudRate
       PathLength
       NumberWavelengths
       NumberTemperatureBins
       TemperatureBins
       cWavelengths
       aWavelengths
       cCleanWaterCalConstant
       aCleanWaterCalConstant
       cTempCompensationVals
       aTempCompensationVals
       aMaxNoise
       cMaxNoise
       aMaxNonConform
       cMaxNonConform
       aMaxDifference
       cMaxDifference
       aMinCounts
       cMinCounts
       rMinCounts
       maxTempSdev
       maxDepthSdev
    end
    properties (SetAccess = private, GetAccess = private)
    end
    methods
     
        %%# ACDeviceFile
         function obj = ACDeviceFile(fileNameIn, fileTypeIn)
            
            % ACDeviceFile creates an object holding all info from a
            %device file
            %
            % SYNOPSIS ACDeviceFile(fileNameIn)
            % INPUT fileNameIn: name of the file to use
            % OUTPUT none
             
             
            % get copy of logger
            L = log4m.getLogger();
             
            if nargin > 0
             

                % check there is a file list              
                if ~isempty( fileNameIn ) 
                    obj.FileName = fileNameIn;
                else
                    L.error('ACDeviceFile.ACDeviceFile', 'Need file name');
                end   
                if ~isempty( fileTypeIn ) 
                    obj.FileType = fileTypeIn;
                else
                    L.error('ACDeviceFile.ACDeviceFile', 'Need file name');
                end                   
                
                fid = fopen(fileNameIn);
                
                if fid == -1
                    L.error('ACDeviceFile.ACDeviceFile', 'Device file not successfully opened')
                else
                      
                    %% Read first 9 lines of file
                    if strcmp(obj.FileType, 'WAP')
                        Header = textscan( fid, '%s', 1, 'Delimiter', '\n');
                        Intro = textscan( fid, '%s', 9, 'Delimiter', '\n');
  
                    else
                        Intro = textscan( fid, '%s', 9, 'Delimiter', '\n');
                    end;

                    % set device name
                    obj.DeviceName = strtrim(Intro{1}{1});

                    % set serial number
                    serialNumber = Intro{1}{2};
                    obj.SerialNumber = strtrim(strtok(serialNumber, ';'));

                    % get version number
                    versionNumber = Intro{1}{3};
                    obj.StructureVersionNumber = strtrim(strtok(versionNumber, ';'));

                    % get calibration temps
                    cals = Intro{1}{4};
%                     tempcals = textscan(cals, '"Tcal:%s C, ical: %s C. The offsets were saved to this file on %s."');
%                     obj.TCal = tempcals{1}{1};
%                     obj.ICal = tempcals{2}{1};
%                     obj.CalDate = strtok(tempcals{3}{1}, '.');
                    tempcals = textscan( cals, '"%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s[\n]');
                    obj.TCal = tempcals{1,2};
                    obj.ICal = tempcals{1,5};
                    obj.CalDate = strtok(tempcals{1,15}, '.');

                    depth = Intro{1}{5};
                    obj.CalibrationDepth = strtok(depth, ';');

                    baudRate = Intro{1}{6};
                    obj.BaudRate = strtrim(strtok(baudRate, ';'));

                    pathLength = Intro{1}{7};
                    obj.PathLength = strtrim(strtok(pathLength, ';'));

                    numWL = Intro{1}{8};
                    obj.NumberWavelengths = strtrim(strtok(numWL, ';'));

                    numTempBins = Intro{1}{9};
                    obj.NumberTemperatureBins = strtrim(strtok(numTempBins, ';'));

                    %%  Read Temperature Bins
                    TempBins = textscan( fid, strcat( repmat('%f', 1, str2num(obj.NumberTemperatureBins)), ' ; temperature bins'), 1, 'Delimiter', '\n\');
                    obj.TemperatureBins = cell2mat(TempBins(:));

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
                    
                    numBins = str2num(obj.NumberTemperatureBins);
                    startNextBin = numBins+1;
                    numTotalBins = (numBins)*2;
                    numRows = str2num(obj.NumberWavelengths);
                    
                    firstcolsSpec = '%s %s %s %10.6f %10.6f';
                    endrowtext = ' "; C and A offset, and C and A temperature correction info"';
                    FormatString = strcat(firstcolsSpec, repmat('%f ', 1, numTotalBins), endrowtext);
                    
                    %% THINK ABOUT CHANGING THIS TO FSCANF IT WILL GO RIGHT INTO MATRIX
                    tempCompValueMatrix = textscan( fid, FormatString, numRows );

                    obj.cWavelengths = cell2mat(tempCompValueMatrix{:,1});
                    obj.aWavelengths = cell2mat(tempCompValueMatrix{:,2});
                    obj.cWavelengths = str2num(obj.cWavelengths(:,2:end));
                    obj.aWavelengths = str2num(obj.aWavelengths(:,2:end));
                    obj.cCleanWaterCalConstant = tempCompValueMatrix{:,4};
                    obj.aCleanWaterCalConstant = tempCompValueMatrix{:,5};
                    tempCompVals = tempCompValueMatrix(:,6:end);
                    
                    
                    
                    obj.cTempCompensationVals = cell2mat(tempCompVals(:, 1:numBins));
                    obj.aTempCompensationVals = cell2mat(tempCompVals(:, startNextBin:end));

                    %% Read Last line of file -- doesn't seem to contain any data
                    FooterFormatString = strcat( repmat('%f ', 1, 11), ' ; ', repmat(' %s ', 1, 11));
                    Footer = textscan( fid, FooterFormatString, 1);
                    obj.aMaxNoise = Footer{:,1};
                    obj.cMaxNoise = Footer{:,2}; 
                    obj.aMaxNonConform = Footer{:,3};
                    obj.cMaxNonConform = Footer{:,4};
                    obj.aMaxDifference = Footer{:,5};
                    obj.cMaxDifference = Footer{:,6};
                    obj.aMinCounts = Footer{:,7};
                    obj.cMinCounts = Footer{:,8};
                    obj.rMinCounts = Footer{:,9};
                    obj.maxTempSdev = Footer{:,10};
                    obj.maxDepthSdev = Footer{:,11};
                    
                    %% Close file
                    closeresult = fclose(fid);
                    if closeresult == 0 
                        L.info('ACDeviceFile.ACDeviceFile', 'File close successful')
                    else
                        L.error('ACDeviceFile.ACDeviceFile', 'File close not successful')
                    end
                
                end   % end if file not opened
 
            else
                L.error('ACDeviceFile.ACDeviceFile', 'Supply an input argument')
            end   % end if nargin > 0
            
        end   % #ACDeviceFile
        
        
        %%# getInfo
        function getInfo(obj)
            
            %# getInfo displays all the information saved in this object
            %#
            %# SYNOPSIS getInfo(obj)
            %# INPUT obj: the object
            %# OUTPUT none
           L = log4m.getLogger();
           L.info('ACDeviceFile.ACDeviceFile', 'FileName:')
           L.info('ACDeviceFile.ACDeviceFile', (obj.FileName))
           L.info('ACDeviceFile.ACDeviceFile', ('Device Name:'))
           L.info('ACDeviceFile.ACDeviceFile', (obj.DeviceName))
           L.info('ACDeviceFile.ACDeviceFile', ('Serial number:'))
           L.info('ACDeviceFile.ACDeviceFile', (obj.SerialNumber))
           
           L.info('ACDeviceFile.ACDeviceFile', ('Structure Version Number:'))
           L.info('ACDeviceFile.ACDeviceFile', (obj.StructureVersionNumber))
           
           L.info('ACDeviceFile.ACDeviceFile', ('Tcal:'))
           L.info('ACDeviceFile.ACDeviceFile', (obj.TCal))
           
           L.info('ACDeviceFile.ACDeviceFile', ('Ical:'))
           L.info('ACDeviceFile.ACDeviceFile', (obj.ICal))
           
           L.info('ACDeviceFile.ACDeviceFile', ('CalDate:'))
           L.info('ACDeviceFile.ACDeviceFile', (obj.CalDate))
           
           L.info('ACDeviceFile.ACDeviceFile', ('CalDepth:'))
           L.info('ACDeviceFile.ACDeviceFile', (obj.CalibrationDepth))
           
           L.info('ACDeviceFile.ACDeviceFile', ('BaudRate:'))
           L.info('ACDeviceFile.ACDeviceFile', (obj.BaudRate))
           
           L.info('ACDeviceFile.ACDeviceFile', ('PathLength:'))
           L.info('ACDeviceFile.ACDeviceFile', (obj.PathLength))
           
           L.info('ACDeviceFile.ACDeviceFile', ('NumberWL:'))
           L.info('ACDeviceFile.ACDeviceFile', (obj.NumberWavelengths))
           
           L.info('ACDeviceFile.ACDeviceFile', ('Number Temp Bins:'))
           L.info('ACDeviceFile.ACDeviceFile', (obj.NumberTemperatureBins))
        
           L.info('ACDeviceFile.ACDeviceFile', ('c wavelengths:'))
           L.info('ACDeviceFile.ACDeviceFile', (obj.cWavelengths))
           
           L.info('ACDeviceFile.ACDeviceFile', ('a wavelengths:'))
           L.info('ACDeviceFile.ACDeviceFile', (obj.aWavelengths))
           
        end   % # getInfo
    end   % # methods
end