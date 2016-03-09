%mapACS_HPLC - Create a GoogleEarth file of the lat/lon coordinates of a 
% cruise, from more than one data set (i.e. from AC-S and HPLC data
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
% Feb 2016; Last revision: 9-Mar-16

%------------- BEGIN CODE --------------
% This is for the ac-s files:

close all;
clear all;

files = dir('c:/users/wendy/documents/data/tara/taramed/processed/final/Tara_ACS_apcp*ap.txt');
dir_acs = 'c:/users/wendy/documents/data/tara/taramed/processed/final';

filelist = {files.name}';
iconStrBase = 'http://maps.google.com/mapfiles/kml/shapes/';
all_acs_tmp =[];

iPrintDateAsName = 1;

for iFiles = 1:length(filelist)

    % create variables
    acs_tmp = [];
    seabassFileName = fullfile(dir_acs, filelist{iFiles});
    shortFileName = filelist{iFiles};
    variablesFromFileName = textscan(shortFileName, 'Tara_ACS_apcp%4s_%3sap_slade.txt');
    year = char(variablesFromFileName{1,1});
    yearday = char(variablesFromFileName{1,2});
    apURL = sprintf('http://misclab.umeoce.maine.edu/images/taramed/%s_%s_wlap.jpg',...
        year, yearday);
    cpURL = sprintf('http://misclab.umeoce.maine.edu/images/taramed/%s_%s_wlcp.jpg',...
        year, yearday);
    ap_html_text = sprintf('ap: <a href="%s"><img src="http://misclab.umeoce.maine.edu/images/taramed/%s_%s_wlap.jpg" height=75 width=75></a>',apURL,year, yearday);
    cp_html_text = sprintf('cp: <a href="%s"><img src="http://misclab.umeoce.maine.edu/images/taramed/%s_%s_wlcp.jpg" height=75 width=75></a>',cpURL,year, yearday);

    datamatrix = importSeaBASS2(seabassFileName);

    date = char(cellstr(datamatrix(:,1)));
    time = char(cellstr(datamatrix(:,2)));
    lat = cell2mat(datamatrix(:,3));
    lon = cell2mat(datamatrix(:,4));
    timestamp = datenum(strcat(date,time),'yyyymmddHH:MM:SS');
    timestamp = datestr(timestamp);

    num_measurements = length(lat);
    
    for iMeasurements=1:num_measurements;
        % for the first measurement of each yearday, designate the sailboat
        % icon
        if iMeasurements == 1
            if iPrintDateAsName == 1
             tmp = ge_point(lon(iMeasurements), lat(iMeasurements), 0,...
                           'iconURL',[iconStrBase,'/sailing.png'],...
                            'msgToScreen',true,...
                            'name', sprintf('%s', timestamp(iMeasurements,1:11)),...
                            'description',...
                        sprintf('<h2>Tara Mediterranean</h2>Yearday: %s<BR> Date:  %s<BR> Data Plots (click to expand):<br> %s <br>%s', yearday, timestamp(iMeasurements,:),ap_html_text, cp_html_text));     
            else
                tmp = ge_point(lon(iMeasurements), lat(iMeasurements), 0,...
                           'iconURL',[iconStrBase,'/sailing.png'],...
                            'msgToScreen',true,...
                            'description',...
                        sprintf('<h2>Tara Mediterranean</h2>Yearday: %s<BR> Date:  %s<BR> Data Plots (click to expand):<br> %s <br>%s', yearday, timestamp(iMeasurements,:),ap_html_text, cp_html_text));     
            end;
            % iterate print name
            if iPrintDateAsName <7
                iPrintDateAsName = iPrintDateAsName + 1;
            else
                iPrintDateAsName = 1;
            end;
        else
            tmp = ge_point(lon(iMeasurements), lat(iMeasurements), 0,...
                           'iconURL',[iconStrBase,'/placemark_circle.png'],...
                           'iconScale', .25,...
                            'msgToScreen',true,...
                            'description',...
                        sprintf('<h2>Tara Mediterranean</h1>Yearday: %s<BR> Date:  %s<BR> Data Plots (click to expand):<br> %s <br>%s', yearday, timestamp(iMeasurements,:),ap_html_text, cp_html_text));                 
        end;
        acs_tmp(end+1:end+length(tmp)) = tmp;
    end;

    all_acs_tmp(end+1:end+length(acs_tmp)) = acs_tmp;

end
    
%% start HPLC data

files_h = dir('c:/users/wendy/documents/data/tara/taramed/HPLC/Tarapigments*csv');
dir_h = 'c:/users/wendy/documents/data/tara/taramed/HPLC';
all_hplc_tmp =[];
filelist = {files_h.name}';
iPrintDateAsName = 1;

for iFiles = 1:length(filelist)
% iFiles=1
    % create variables
    hplc_tmp = [];
    seabassFileName = fullfile(dir_h, filelist{iFiles});
    shortFileName = filelist{iFiles};
%     variablesFromFileName = textscan(shortFileName, 'Tarapigments_2014_%s.csv');

    datamatrix = importHPLC(seabassFileName);

    date = char(cellstr(datamatrix(:,1)));
    time = char(cellstr(datamatrix(:,2)));
    lat = str2double(datamatrix(:,3));
    lon = str2double(datamatrix(:,4));
    timestamp = datenum(strcat(date,time),'yyyymmddHH:MM:SS');
    timestamp = datestr(timestamp);

    num_measurements = length(lat);
    
    for iMeasurements=1:num_measurements;
        tmp = ge_point(lon(iMeasurements), lat(iMeasurements), 0,...
           'iconURL',[iconStrBase,'/placemark_circle_highlight.png'],...
            'msgToScreen',true,...
            'name', sprintf('%s', timestamp(iMeasurements,1:11)),...
            'description',...
            sprintf('<h2>Tara Mediterranean HPLC</h2>Date:  %s', timestamp(iMeasurements,:)));     

        hplc_tmp(end+1:end+length(tmp)) = tmp;
    end;

    all_hplc_tmp(end+1:end+length(hplc_tmp)) = hplc_tmp;

end % for iFiles

%% Merge both together
all_tmp = [all_acs_tmp all_hplc_tmp];
% ---------------------------------------------------------------------------------------
kmlStr01 = ge_folder('ACS Data', [all_tmp]);
ge_output('C:\Users\Wendy\Dropbox\Tara_Med_acshplc.kml',kmlStr01);

%------------- END OF CODE --------------