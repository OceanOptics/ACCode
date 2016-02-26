%mapACS - Create a GoogleEarth file of the lat/lon coordinates of a cruise
%
% Other m-files required: Google Earth Toolbox
% Subfunctions: importSeaBASS2
% MAT-files required: none
%
% See also: ge_point, importSeaBASS2

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% Feb 2016; Last revision: 15-FEb-16

%------------- BEGIN CODE --------------

close all;
clear all;

files = dir('c:/users/wendy/documents/data/tara/taraarctic/ag/*.csv');
dir = 'c:/users/wendy/documents/data/tara/taraarctic/ag';

filelist = {files.name}';
iconStrBase = 'http://maps.google.com/mapfiles/kml/shapes/';
all_acs_tmp =[];

for iFiles = 1:length(filelist)

    % create variables
    acs_tmp = [];
    seabassFileName = fullfile(dir, filelist{iFiles});
    shortFileName = filelist{iFiles};
    variablesFromFileName = textscan(shortFileName, 'Tara_ACS_apcp%4s%s');
    year = char(variablesFromFileName{1,1});
%     yearday = char(variablesFromFileName{1,2});

    datamatrix = importfileSeaBASSArctic(seabassFileName);

    date = char(cellstr(datamatrix(:,1)));
    time = char(cellstr(datamatrix(:,2)));
    lat = str2double(datamatrix(:,3));
    lon = str2double(datamatrix(:,4));
    timestamp = datenum(strcat(date,time),'yyyymmddHH:MM:SS');
    timestamp = datestr(timestamp);

    num_measurements = length(lat);
    
    for iMeasurements=1:num_measurements;
        % for the first measurement of each yearday, designate the sailboat
        % icon
        if iMeasurements == 1
             tmp = ge_point(lon(iMeasurements), lat(iMeasurements), 0,...
                           'iconURL',[iconStrBase,'/sailing.png'],...
                            'msgToScreen',true,...
                            'description',...
                        sprintf('<h2>Tara Mediterranean</h2>Date:  %s<BR>', timestamp(iMeasurements,:)));               
        else
            tmp = ge_point(lon(iMeasurements), lat(iMeasurements), 0,...
                           'iconURL',[iconStrBase,'/placemark_circle.png'],...
                           'iconScale', .25,...
                            'msgToScreen',true,...
                            'description',...
                        sprintf('<h2>Tara Mediterranean</h2>Date:  %s<BR>', timestamp(iMeasurements,:)));  
        end;
        acs_tmp(end+1:end+length(tmp)) = tmp;
    end;

    all_acs_tmp(end+1:end+length(acs_tmp)) = acs_tmp;

end
    
kmlStr01 = ge_folder('ACS Data',[all_acs_tmp]);
ge_output('C:\Users\Wendy\Dropbox\Tara_Arctic_acs.kml',kmlStr01);

%------------- END OF CODE --------------