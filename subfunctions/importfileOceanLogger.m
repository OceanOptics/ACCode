function [timestamp,airtemp,humidity,par1,tir1,airtemp2,humidity2,par2,tir2,baro1,baro2,tstemp,conductivity,salinity,sound_velocity,chlorophyll,sampletemp,flowrate,sstemp,trans,sstemp2] = importfileOceanLogger(filename, startRow, endRow)
%IMPORTFILE1 Import numeric data from a text file as column vectors.
%   [TIMESTAMP,AIRTEMP,HUMIDITY,PAR1,TIR1,AIRTEMP2,HUMIDITY2,PAR2,TIR2,BARO1,BARO2,TSTEMP,CONDUCTIVITY,SALINITY,SOUND_VELOCITY,CHLOROPHYLL,SAMPLETEMP,FLOWRATE,SSTEMP,TRANS,SSTEMP2]
%   = IMPORTFILE1(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [TIMESTAMP,AIRTEMP,HUMIDITY,PAR1,TIR1,AIRTEMP2,HUMIDITY2,PAR2,TIR2,BARO1,BARO2,TSTEMP,CONDUCTIVITY,SALINITY,SOUND_VELOCITY,CHLOROPHYLL,SAMPLETEMP,FLOWRATE,SSTEMP,TRANS,SSTEMP2]
%   = IMPORTFILE1(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [timestamp,airtemp,humidity,par1,tir1,airtemp2,humidity2,par2,tir2,baro1,baro2,tstemp,conductivity,salinity,sound_velocity,chlorophyll,sampletemp,flowrate,sstemp,trans,sstemp2] = importfile1('oceanlogger.298',1, 17205);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/09/13 12:23:59

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format string for each line of text:
%   column5: text (%s)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%*s%*s%*s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
timestamp = dataArray{:, 1};
airtemp = dataArray{:, 2};
humidity = dataArray{:, 3};
par1 = dataArray{:, 4};
tir1 = dataArray{:, 5};
airtemp2 = dataArray{:, 6};
humidity2 = dataArray{:, 7};
par2 = dataArray{:, 8};
tir2 = dataArray{:, 9};
baro1 = dataArray{:, 10};
baro2 = dataArray{:, 11};
tstemp = dataArray{:, 12};
conductivity = dataArray{:, 13};
salinity = dataArray{:, 14};
sound_velocity = dataArray{:, 15};
chlorophyll = dataArray{:, 16};
sampletemp = dataArray{:, 17};
flowrate = dataArray{:, 18};
sstemp = dataArray{:, 19};
trans = dataArray{:, 20};
sstemp2 = dataArray{:, 21};

