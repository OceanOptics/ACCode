%% OutputManager - SCRIPT which takes processed data, saved as an 
% ProcessedData object called "pd" and outputs SeaBASS data files and plots
%
%
% Requires: the .mat file saved in ProcessingManager, if not in memory
%
% Outputs:  SeaBASS txt files
%
% Other m-files required: readIngestParameters, ProcessedData, 
% Subfunctions: checkTimestamps
% MAT-files required: "acsPROC_(YEAR)_(DAY)" saved by PreProcessingManager,
%  "acsPREPROC_(YEAR)_(DAY)" saved by PreProcessingManager
%
% See also: readIngestParameters, ProcessedData, checkTimestamps
% 
% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% May 2015; Last revision: 16-Feb-16
% 

%----------------------------- BEGIN CODE ---------------------------------
%%Load data file from disk if necessary

if params.RUN.LOAD_OUTPUT_DATA_FROM_DISK
    
    % load ap and cp data
    matFileName = fullfile(params.INGEST.DATA_INPUT_DIRECTORY, ...
        strcat('acsPROC', '_', num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY)));
    load(matFileName);

    % load ancillary data
    matFileName = fullfile(params.INGEST.DATA_INPUT_DIRECTORY, ...
        strcat('acsPREPROC', '_', num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY)));
    load(matFileName);  
    
end;

   
%% ----------------------------------------------------------------------- 
% SECTION 1: OUTPUT SEABASS FILES
% -----------------------------------------------------------------------

% ----------------------------------------------------------------------
% PART 1: GET DATA
% ----------------------------------------------------------------------

% get variables needed for output:
% ap, cp, ap_std, cp_std
% get timestamps
ap_timestamps = pd.getVar('name', 'ap', 'data', 'timestamps');
cp_timestamps = pd.getVar('name', 'cp', 'data', 'timestamps');

% get ap - for each correction - more than one correction is possible
if params.PROCESS.SCATTERING_CORR_SLADE
%     ap_data_slade = pd.getVar('name', 'ap', 'data', 'data_slade', 'level', 'corrected');
    ap_data_slade = pd.getVar('name', 'ap', 'data', 'data_slade');
end;
if params.PROCESS.SCATTERING_CORR_ROTTGERS
    ap_data_rottgers = pd.getVar('name', 'ap', 'data', 'data_rottgers', 'level', 'corrected');
end;
if params.PROCESS.SCATTERING_CORR_FLAT
    ap_data_flat = pd.getVar('name', 'ap', 'data', 'data_flat', 'level', 'corrected');
end;

% get cp
% cp_data = pd.getVar('name', 'cp', 'data', 'data', 'level', 'corrected');
cp_data = pd.getVar('name', 'cp', 'data', 'data');

% get ap uncertainty
% there can be more than one correction calculated - but we can only report
% one to seabass
% THIS IS THE SECTION THAT NEEDS TO BE CHANGED
% if params.PROCESS.AP_UNCERTAINTY_BETWEEN_CORRECTIONS
if params.OUTPUT.USE_AP_UNCERTAINTY_BETWEEN_CORRECTIONS
    ap_uncertainty = pd.getVar('name', 'ap', 'data', 'uncertainty_between_corrections', 'level', 'corrected');
elseif params.OUTPUT.USE_AP_UNCERTAINTY_1
    ap_uncertainty = pd.getVar('name', 'ap', 'data', 'uncertainty', 'level', 'corrected');
elseif params.OUTPUT.USE_AP_STD
    % NEED STATEMENT HERE ONCE OTHER CODE CHANGES
end;
    
% get cp uncertainty
if params.OUTPUT.USE_AP_UNCERTAINTY_1
    cp_uncertainty = pd.getVar('name', 'cp', 'data', 'uncertainty');
elseif params.OUTPUT.USE_AP_STD
    % NEED STATEMENT HERE ONCE OTHER CODE CHANGES
end;


% --------------------------------------------------------------
% added 2/25/16 for TaraMed Hack:
ap_std = pd.getVar('name','aTSW','data','std','level','binned');
cp_std = pd.getVar('name','cTSW','data','std','level','binned');

% create an index of timestamps to remove
% (This was done to other data previously)
timestamps_a = pd.var.aTSW.L3.binnedTime;
timestamps_c = pd.var.cTSW.L3.binnedTime;
[r,c] = size(timestamps_a);
if r<c
    timestamps_a = timestamps_a';
end;
[r,c] = size(timestamps_c);
if r<c
    timestamps_c = timestamps_c';
end;
% check ap & cp timestamps are the same, otherwise make the same
a_ts_size = size(timestamps_a);
c_ts_size = size(timestamps_c);
if a_ts_size(1) ~= c_ts_size(1)
    L.info('ProcessingManager','ap size and cp size diff')
    if abs(c_ts_size(1) - a_ts_size(1)) == 1
        % only one timestamp in differnece - difference of seconds
        if c_ts_size(1) > a_ts_size(1)
            %cp is bigger
            % if end-1 for cp is same as end for ap
            if sum(datevec(timestamps_a(end)) == datevec(timestamps_c(end-1))) == 6
                timestamps_c(end,:) = [];
                cp_std(end,:) = [];
            else
                L.error('ProcessingManager','gap between a tiimestamps and c timestamps is larger than 1');
            end;
        elseif a_ts_size(1) > c_ts_size(1)
            %ap is bigger
            % if end-1 for ap is same as end for cp
            if sum(datevec(timestamps_a(end-1)) == datevec(timestamps_c(end))) == 6
                timestamps_a(end) = [];
                ap_std(end) = [];
            else
                L.error('ProcessingManager','gap between a tiimestamps and c timestamps is larger than 1');
            end;
        else
            L.error('ProcessingManager','neither a nor c seems to be bigger');
        end;
    else
        L.error('ProcessingManager','gap between a tiimestamps and c timestamps is larger than 1');
    end;
end;
% check again
ap_size2 = size(pd.var.ap.L5.timestamps);
cp_size2 = size(pd.var.cp.L5.timestamps);
if ap_size2(1) ~= cp_size2(1)
    L.error('ProcessingManager', 'tried to fix size gap between a and c tiemstamps and failed')
else
    L.info('ProcessingManager', 'adjusted gap between and c timestamps by 1');
end;
% ---------------------------------------


removeIndex = timestamps_a > datenum(pd.meta.Params.INGEST.YEAR, ...
    pd.meta.Params.INGEST.MONTH, pd.meta.Params.INGEST.DAY, 23, 59, 59) | ...
    timestamps_a < datenum(pd.meta.Params.INGEST.YEAR, ...
    pd.meta.Params.INGEST.MONTH, pd.meta.Params.INGEST.DAY);
timestamps_a_after = removerows(timestamps_a, removeIndex);
ap_std = removerows(ap_std, removeIndex);
% ap_std(removeIndex,:) = [];



removeIndex_c = timestamps_c > datenum(pd.meta.Params.INGEST.YEAR, ...
    pd.meta.Params.INGEST.MONTH, pd.meta.Params.INGEST.DAY, 23, 59, 59) | ...
    timestamps_c < datenum(pd.meta.Params.INGEST.YEAR, ...
    pd.meta.Params.INGEST.MONTH, pd.meta.Params.INGEST.DAY);
timestamps_c_after = removerows(timestamps_c, removeIndex_c);
cp_std = removerows(cp_std, removeIndex_c);
% --------------------------------------------------------------
% added 2/25/16 for TaraMed Hack:
%Hack c wavelenght cutoff at 750
%             wavelengths = obj.meta.DeviceFile.(thisWL);
%             data = obj.getVar('name', thisType,'data','data');
c_orig_wl = pd.meta.DeviceFile.cWavelengths;
c_cut_wl = pd.var.cp.L6.wavelengths;

cp_std_EndAt750 = interp1( c_orig_wl, cp_std', c_cut_wl, 'linear', 'extrap');
cp_std_EndAt750 = cp_std_EndAt750';

% --------------------------------------------------------------
% added 2/25/16 for TaraMed Hack:
% hack interp a onto c
a_orig_wl = pd.meta.DeviceFile.aWavelengths;
a_cut_wl = pd.var.ap.L7.wavelengths;

ap_std_MatchedUnc = interp1( a_orig_wl, ap_std', c_cut_wl, 'linear', 'extrap');
ap_std_MatchedUnc = ap_std_MatchedUnc';

% copy into final vars
ap_std_to_print = ap_std_MatchedUnc;
cp_std_to_print = cp_std_EndAt750;

% --------------------------------------------------------------
% --------------------------------------------------------------
% GET ANCILLARY DATA: Temperature, Salinity, GPS
temp_data = allData.TemperatureData.var.L3.BinnedLabTempData;
temp_timestamps = allData.TemperatureData.var.L3.BinnedTimestamps;
sal_timestamps = allData.SalinityData.var.L3.BinnedTimestamps;
sal_data = allData.SalinityData.var.L3.BinnedData;
lat_data = allData.GPSData.var.L3.BinnedLatData;
lon_data = allData.GPSData.var.L3.BinnedLonData;
gps_timestamps = allData.GPSData.var.L3.BinnedTimestamps;

% should all be the same
wavelengths = pd.var.ap.L8.wavelengths_slade;


%% ----------------------------------------------------------------------
% PART 2: CHECK/MANIPULATE DATA
% ----------------------------------------------------------------------
%
% check timestamps are all the same
% CHANGED 10/22
if checkTimestamps(ap_timestamps, cp_timestamps)
    L.debug('OutputManager', 'a and c timestamps same')
else
    L.error('OutputManager', 'a and c timestamps diff')
end;
if checkTimestamps(ap_timestamps, temp_timestamps)
    L.debug('OutputManager', 'a and temp timestamps same')
else
    L.error('OutputManager', 'a and temp timestamps diff')
end;
if checkTimestamps(ap_timestamps, gps_timestamps)
    L.debug('OutputManager', 'a and gps timestamps same')
else
    L.error('OutputManager', 'a and gps timestamps diff')
end;
if checkTimestamps(ap_timestamps, sal_timestamps)
    L.debug('OutputManager', 'a and sal timestamps same')
else
    L.error('OutputManager', 'a and sal timestamps diff')
end;

%%
goodrows = logical(isfinite(ap_data_slade(:,20)) .* isfinite(lat_data) .* ...
    isfinite(cp_data(:,20)));

%%
% Replace values below 0.005 with -9999
replaceWith999Index = cp_data < -0.005;
numReplaced = sum(replaceWith999Index);
totalNumReplaced = sum(numReplaced);
L.info('OutputManager', sprintf('Replaced %u values in cp_data with -9999', totalNumReplaced)); 

cp_data(replaceWith999Index) = -9999;
% cp_uncertainty(replaceWith999Index) = -9999;
cp_std_to_print(replaceWith999Index) = -9999;

if params.PROCESS.SCATTERING_CORR_SLADE
    replaceWith999Index = ap_data_slade < -0.005;
    numReplaced = sum(replaceWith999Index);
    totalNumReplaced = sum(numReplaced);
    L.info('OutputManager', sprintf('Replaced %u values in ap_data_slade with -9999', totalNumReplaced));
    ap_data_slade(replaceWith999Index) = -9999;
%     ap_uncertainty(replaceWith999Index) = -9999;
    ap_std_to_print(replaceWith999Index) = -9999;
end

% if params.PROCESS.SCATTERING_CORR_ROTTGERS
%     replaceWith999Index = ap_data_rottgers < -0.005;
%     numReplaced = sum(replaceWith999Index);
%     totalNumReplaced = sum(numReplaced);
%     L.info('OutputManager', sprintf('Replaced %u values in ap_data_rottgers with -9999', totalNumReplaced));
%     ap_data_rottgers(replaceWith999Index) = -9999;
%     ap_uncertainty(replaceWith999Index) = -9999;
% end

% if params.PROCESS.SCATTERING_CORR_FLAT
%     replaceWith999Index = ap_data_flat < -0.005;
%     numReplaced = sum(replaceWith999Index);
%     totalNumReplaced = sum(numReplaced);
%     L.info('OutputManager', sprintf('Replaced %u values in ap_data_flat with -9999', totalNumReplaced));    
%     ap_data_flat(replaceWith999Index) = -9999;
%     ap_uncertainty(replaceWith999Index) = -9999;
% end

% bin_tempc(isnan(bin_tempc)) = -9999;
replaceWith999Index = isnan(temp_data);
numReplaced = sum(replaceWith999Index);
totalNumReplaced = sum(numReplaced);
L.info('OutputManager', sprintf('Replaced %u values in temp_data with -9999', totalNumReplaced));  
temp_data(replaceWith999Index) = -9999;

replaceWith999Index = isnan(sal_data);
numReplaced = sum(replaceWith999Index);
totalNumReplaced = sum(numReplaced);
L.info('OutputManager', sprintf('Replaced %u values in sal_data with -9999', totalNumReplaced));  
sal_data(replaceWith999Index) = -9999;

% get start and end dates and times
% same for all files for this yearday
timestamps = cellstr(datestr(ap_timestamps));
temp_date = datestr(datenum(timestamps),'yyyymmdd');
temp_time = datestr(datenum(timestamps),'HH:MM:SS');
timestamps_date_to_print = cellstr(temp_date(goodrows,:));
timestamps_time_to_print = cellstr(temp_time(goodrows,:));

% same for all files for this yearday - NEEDS TO BE DONE AFTER BADROWS
% FILTERED OUT
% start_date = datestr(datenum(timestamps(1, 1)),'yyyymmdd');
% end_date = datestr(datenum(timestamps(end, 1)),'yyyymmdd');
% start_time = datestr(datenum(timestamps(1, 1)),'HH:MM:SS');
% end_time = datestr(datenum(timestamps(end, 1)),'HH:MM:SS');
start_date = timestamps_date_to_print{1,1};
end_date = timestamps_date_to_print{end,1};
start_time = timestamps_time_to_print{1,1};
end_time = timestamps_time_to_print{end,1};

% same for all files for this yearday
% northLat = max(lat_data);
% southLat = min(lat_data);
% eastLon = max(lon_data);
% westLon = min(lon_data);

% ----------------------------------------------------------------------
% PART 2: PRODUCE HEADER
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% Step 1: Set up variables
% ----------------------------------------------------------------------
sb_fname_xls = [params.OUTPUT.SEABASS_FILE_PREFIX num2str(params.INGEST.YEAR) '_' num2str(params.INGEST.YEAR_DAY) '.xls'];
sb_fname_ascii = [params.OUTPUT.SEABASS_FILE_PREFIX num2str(params.INGEST.YEAR) '_' num2str(params.INGEST.YEAR_DAY)];
cruise = strrep(params.INGEST.CRUISE_LEG, ' ', '');

% cp_std_cutoff = {ap_data_slade, ap_data_rottgers, cp_data, ap_uncertainty, cp_uncertainty};
% dataFiles = {'ap_slade', 'ap_rottgers', 'cp', 'ap_uncertainty', 'cp_uncertainty'};
% wlPrefix = {'ap', 'ap', 'cp', 'ap', 'cp'};
if params.OUTPUT.USE_SLADE
    data_type = {ap_data_slade, cp_data};
elseif params.OUTPUT.USE_ROTTGERS
    data_type = {ap_data_rottgers, cp_data};
elseif params.OUTPUT.USE_FLAT
    data_type = {ap_data_flat, cp_data};
else
    L.error('OutputManager','No ap data chosen for output');
end;

std_type = {ap_std_to_print, cp_std_to_print};
dataFiles = {'ap', 'cp'};
wlPrefix = {'ap', 'cp'};

if ~exist(params.INGEST.DATA_OUTPUT_DIRECTORY, 'dir')
    mkdir(params.INGEST.DATA_OUTPUT_DIRECTORY)
    L.info('OutputManager','made new direcotry for output');
end;
%%
for iData = 1:length(data_type)
    % ----------------------------------------------------------------------
    % Step 2: Create Header
    % ----------------------------------------------------------------------
    thisData = data_type{iData}; %ap_data_slade;
    thisSTD = std_type{iData};
    thisFileName = dataFiles{iData}; %;
    thisPrefix = wlPrefix{iData}; %;

    data_to_print = [lat_data lon_data temp_data sal_data thisData thisSTD];
    data_to_print = data_to_print(goodrows,:);

    northLat = max(data_to_print(:,1));
    southLat = min(data_to_print(:,1));
    eastLon = max(data_to_print(:,2));
    westLon = min(data_to_print(:,2));
    
    % create list of units, get rid of last ','
    apcp_units = char(repmat(uint8(strcat(params.INGEST.AC_UNITS, ',')), 1, 2*(length(wavelengths))));
    apcp_units(length(apcp_units)) = '';

    % build the extension name from the dataFile
    extension = strcat(thisFileName, '.txt');
    seabassFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, strcat(sb_fname_ascii, extension));

    sb_hdr = {'date', 'time', 'lat', 'lon', 'Wt', 'sal'};
    % put all these fields into a matrix
    sb_hdr = [sb_hdr strcat(thisPrefix, cellstr(num2str(wavelengths)))' strcat(thisPrefix, cellstr(num2str(wavelengths)), '_sd')'];
    % remove whitespace
    sb_hdr = strrep(sb_hdr, ' ', '');

    % name file
    % fpath='d:\misclab\tara\Tara Processing\working\';
%     fid = fopen(seabassFileName, 'wt');
    fid = fopen(seabassFileName, 'w', 'n', 'US-ASCII');
    
    fprintf(fid,'/begin_header\n');
    fprintf(fid,'/investigators=%s\n', params.OUTPUT.INVESTIGATOR );
    fprintf(fid,'/affiliations=%s\n', params.OUTPUT.AFFILIATION );
    fprintf(fid,'/contact=%s\n', params.OUTPUT.CONTACT);
    fprintf(fid,'/experiment=%s\n', params.OUTPUT.EXPERIMENT);
    fprintf(fid,'/cruise=');
    fprintf(fid,cruise);
    fprintf(fid,'\n');
    fprintf(fid,'/station=NA\n');
    fprintf(fid,'/data_file_name=');
    fprintf(fid,strcat(sb_fname_ascii,extension));
    fprintf(fid,'\n');
    fprintf(fid,'/documents=%s\n', params.OUTPUT.DOCUMENTATION);
    fprintf(fid,'/calibration_files=%s\n', params.OUTPUT.CALIBRATION_FILES);
    fprintf(fid,'/data_type=%s\n', params.OUTPUT.DATA_TYPE);
    fprintf(fid,'/data_status=%s\n', params.OUTPUT.DATA_STATUS);
    fprintf(fid, '/start_date=%s\n', start_date);
    fprintf(fid, '/end_date=%s\n', end_date);
    fprintf(fid,'/start_time=%s[GMT]\n', start_time);
    fprintf(fid,'/end_time=%s[GMT]\n', end_time);
    fprintf(fid,'/north_latitude=%5.3f[DEG]\n', northLat);
    fprintf(fid,'/south_latitude=%5.3f[DEG]\n', southLat);
    fprintf(fid,'/east_longitude=%5.3f[DEG]\n', eastLon);
    fprintf(fid,'/west_longitude=%5.3f[DEG]\n', westLon);
    fprintf(fid,'/water_depth=NA\n');
    fprintf(fid,'/measurement_depth=%s\n', params.OUTPUT.MEASUREMENT_DEPTH); % not allowed if depth is in the data
    fprintf(fid,'/secchi_depth=NA\n');
    fprintf(fid,'/cloud_percent=NA\n');
    fprintf(fid,'/wind_speed=NA\n');
    fprintf(fid,'/wave_height=NA\n');
    fprintf(fid,'/missing=-9999\n');
    fprintf(fid,'/delimiter=space\n');

    fprintf(fid,'/fields=');
    for j=1:length(sb_hdr)
        if j < length(sb_hdr)
            fprintf(fid, '%s,',cell2mat(sb_hdr(j)));
            % added:
        else
            fprintf(fid, '%s',cell2mat(sb_hdr(j)));
        end;
    end
    fprintf(fid,'\n');
%        fprintf(fid,['/units=yyyymmdd,hh:mm:ss,degrees,degrees,degreesC,PSU,' apcp_units '\n']);

    fprintf(fid,['/units=yyyymmdd,hh:mm:ss,degrees,degrees,degreesC,PSU,' apcp_units '\n']);
    fprintf(fid,'/end_header\n');

    % ----------------------------------------------------------------------
    % Step 3: Create Data
    % ----------------------------------------------------------------------

    % for i=1:m
    %     fprintf(fid, '%8s ',cell2mat(tmp_date(i,1)));
    %     fprintf(fid, '%8s ',cell2mat(tmp_time(i,1)));
    %     % fprintf(fid, '%20s ',cell2mat(sb_datetime(i)));
    %     fprintf(fid,'%6.4f ',sb_dat_ap(i,:));
    %     fprintf(fid,'\n');
    % end
    [nRows, ~]=size(data_to_print);
    for thisRow = 1:nRows
        fprintf(fid, '%8s ', cell2mat( timestamps_date_to_print( thisRow,1 )));
        fprintf(fid, '%8s ', cell2mat( timestamps_time_to_print( thisRow,1 )));
        % fprintf(fid, '%20s ',cell2mat(sb_datetime(i)));
%         fprintf(fid,'%6.4f ', data_to_print(thisRow,:));
        % print all columns except last with a space after
        fprintf(fid,'%6.4f ', data_to_print(thisRow,1:end-1));
        % print last column without a space
        fprintf(fid,'%6.4f', data_to_print(thisRow,end));
        fprintf(fid,'\n');
    end
    fclose(fid);
    
end   % for loop for data


%% ---------------------------------------------------------------------- 
% SECTION 3: CREATE FINAL PLOTS
% -----------------------------------------------------------------------
%%
fignum = 41;
wavelengthToPlot = 20;

figure(fignum);
ax1 = subplot(2,1,1);
hold on;
grid on;
plot(pd.var.c.L1.timestamps, pd.var.c.L1.data(:, wavelengthToPlot), 'b');
scatter(pd.var.cFSW.L3.binnedTime, pd.var.cFSW.L3.median(:,20), 5, 'o', 'fill', 'MarkerFaceColor', 'green', 'MarkerEdgeColor','green');
xlabel('Timestamps');
ylabel('c Data');
dynamicDateTicks;
xlim([ datenum(params.INGEST.YEAR, params.INGEST.MONTH, params.INGEST.DAY, 0, 0, 0) ,...
       datenum(params.INGEST.YEAR, params.INGEST.MONTH, params.INGEST.DAY + 1, 0, 0, 0) ] );
% ylim([-0.1, .5])
% legend('Synchronized c data', 'FSW Bins', 'Interpolated data');
legend('Synchronized c data', 'FSW Bins');
title_text = sprintf('Filtered c Data - using median\n Cruise: %s Leg: %s \n%u-%s-%u (Yearday: %u)', ...
    params.INGEST.CRUISE, params.INGEST.CRUISE_LEG, params.INGEST.DAY, ...
    params.INGEST.MONTH_TEXT, params.INGEST.YEAR, params.INGEST.YEAR_DAY);    
title(title_text,'fontsize',12);

ax2 = subplot(2,1,2);
hold on;
grid on;
plot(pd.var.a.L1.timestamps, pd.var.a.L1.data(:, wavelengthToPlot), 'b');
scatter(pd.var.aFSW.L3.binnedTime, pd.var.aFSW.L3.median(:,20), 5, 'o', 'fill', 'MarkerFaceColor', 'green', 'MarkerEdgeColor','green');
xlabel('Timestamps');
ylabel('a Data');
dynamicDateTicks;
xlim([ datenum(params.INGEST.YEAR, params.INGEST.MONTH, params.INGEST.DAY, 0, 0, 0) ,...
       datenum(params.INGEST.YEAR, params.INGEST.MONTH, params.INGEST.DAY + 1, 0, 0, 0) ] );
% ylim([ -0.05, 0.15])
legend('Synchronized a data', 'FSW Bins');
title_text = sprintf('Filtered a Data - using median\n Cruise: %s Leg: %s \n%u-%s-%u (Yearday: %u)', ...
    params.INGEST.CRUISE, params.INGEST.CRUISE_LEG, params.INGEST.DAY, ...
    params.INGEST.MONTH_TEXT, params.INGEST.YEAR, params.INGEST.YEAR_DAY);
title(title_text,'fontsize',12);

linkaxes([ax1, ax2], 'x')
saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, strcat(num2str(params.INGEST.YEAR), ...
    '_', num2str(params.INGEST.YEAR_DAY),  '_filtered')));
saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, strcat(num2str(params.INGEST.YEAR), ...
    '_', num2str(params.INGEST.YEAR_DAY), '_filtered.jpg')));
    
%% FIGURE 2: "PEACHY"
fignum = 42;
figure(fignum);
ax1 = subplot(2,1,1);
hold on;
grid on;
plot(pd.var.c.L1.timestamps, pd.var.c.L1.data(:, wavelengthToPlot), 'k');
scatter(pd.var.cTSW.L3.binnedTime, pd.var.cTSW.L3.median(:,20), 5, 'o', 'MarkerFaceColor', 'green', 'MarkerEdgeColor','green');
scatter(pd.var.cFSW.L4.binnedTime, pd.var.cFSW.L4.interpolatedData(:,20),10,'o','fill','MarkerEdgeColor','blue');
xlabel('Timestamps');
ylabel('c Data');
dynamicDateTicks;
xlim([ datenum(params.INGEST.YEAR, params.INGEST.MONTH, params.INGEST.DAY, 0, 0, 0) ,...
       datenum(params.INGEST.YEAR, params.INGEST.MONTH, params.INGEST.DAY + 1, 0, 0, 0) ] );
% ylim([ -0.05, 0.5])
legend('Synchronized c data', 'TSW Bins', 'Interpolated data');
    title_text = sprintf('Interpolated Filtered c Data - using median\n Cruise: %s Leg: %s \n%u-%s-%u (Yearday: %u)', ...
    params.INGEST.CRUISE, params.INGEST.CRUISE_LEG, params.INGEST.DAY, ...
    params.INGEST.MONTH_TEXT, params.INGEST.YEAR, params.INGEST.YEAR_DAY);
title(title_text,'fontsize',12);

% dynamicDateTicks;
ax2 = subplot(2,1,2);
hold on;
grid on;
plot(pd.var.a.L1.timestamps, pd.var.a.L1.data(:, wavelengthToPlot), 'k');
scatter(pd.var.aTSW.L3.binnedTime, pd.var.aTSW.L3.median(:,20), 5, 'o', 'fill', 'MarkerFaceColor', 'green', 'MarkerEdgeColor','green');
scatter(pd.var.aFSW.L4.binnedTime, pd.var.aFSW.L4.interpolatedData(:,20), 5, 'fill', 'MarkerEdgeColor', 'blue');
xlabel('Timestamps');
ylabel('a Data');
dynamicDateTicks;
xlim([ datenum(params.INGEST.YEAR, params.INGEST.MONTH, params.INGEST.DAY, 0, 0, 0) ,...
       datenum(params.INGEST.YEAR, params.INGEST.MONTH, params.INGEST.DAY + 1, 0, 0, 0) ] );
% ylim([ -0.03, 0.5])
legend('Synchronized a data', 'TSW Bins', 'Interpolated data');
    title_text = sprintf('Interpolated Filtered a Data - using median\n Cruise: %s Leg: %s \n%u-%s-%u (Yearday: %u)', ...
    params.INGEST.CRUISE, params.INGEST.CRUISE_LEG, params.INGEST.DAY, ...
    params.INGEST.MONTH_TEXT, params.INGEST.YEAR, params.INGEST.YEAR_DAY);
title(title_text,'fontsize',12);
% 
 linkaxes([ax1, ax2], 'x')

saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
    strcat(num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY), '_interp')));
saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
    strcat(num2str(params.INGEST.YEAR), '_',num2str(params.INGEST.YEAR_DAY), '_interp.jpg')));    

%% FIGURE 3: wlap
fignum=43;
figure(fignum)
hold on;
grid on;
plot(pd.var.ap.L8.wavelengths_slade, pd.var.ap.L9.data_slade);
xlabel('Wavelength')
ylabel('ap - corrected using Slade')
title_text = sprintf('Corrected ap vs. Wavelength\n Cruise: %s Leg: %s \n%u-%s-%u (Yearday: %u)', ...
    params.INGEST.CRUISE, params.INGEST.CRUISE_LEG, params.INGEST.DAY, ...
    params.INGEST.MONTH_TEXT, params.INGEST.YEAR, params.INGEST.YEAR_DAY);
title(title_text, 'fontsize', 12);

saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
    strcat(num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY), '_wlap')));
saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
    strcat(num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY), '_wlap.jpg')));        

 %% FIGURE 4: wlap
fignum=44;
figure(fignum)
hold on;
grid on;
plot(pd.var.ap.L8.wavelengths_rottgers, pd.var.ap.L9.data_rottgers);
xlabel('Wavelength')
ylabel('ap - corrected using Rottgers')
title_text = sprintf('Corrected ap vs. Wavelength\n Cruise: %s Leg: %s \n%u-%s-%u (Yearday: %u)', ...
    params.INGEST.CRUISE, params.INGEST.CRUISE_LEG, params.INGEST.DAY, ...
    params.INGEST.MONTH_TEXT, params.INGEST.YEAR, params.INGEST.YEAR_DAY);
title(title_text, 'fontsize', 12);

saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
    strcat(num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY), '_wlap_rottgers')));
saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
    strcat(num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY), '_wlap_rottgers.jpg')));     

%% FIGURE 5: wlcp
fignum=45;
figure(fignum)
hold on;
grid on;
plot(pd.var.cp.L8.wavelengths, pd.var.cp.L9.data);
xlabel('Wavelength')
ylabel('cp - corrected')
title_text = sprintf('Corrected cp vs. Wavelength\n Cruise: %s Leg: %s \n%u-%s-%u (Yearday: %u)', ...
    params.INGEST.CRUISE, params.INGEST.CRUISE_LEG, params.INGEST.DAY, ...
    params.INGEST.MONTH_TEXT, params.INGEST.YEAR, params.INGEST.YEAR_DAY);
title(title_text, 'fontsize', 12);  

saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
    strcat(num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY), '_wlcp')));
saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
    strcat(num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY), '_wlcp.jpg'))); 
    

%% write params to .ini file
struct2ini( fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
    strcat(num2str(params.INGEST.YEAR), '_',num2str(params.INGEST.YEAR_DAY), '.ini')), params );
   
%% close variables

if params.INGEST.CLEAR_VARS
    clear ans;
    clear ap_data_rottgers;
    clear ap_data_slade;
    clear ap_timestamps;
    clear ap_uncertainty;
    clear apcp_units;
    clear ax1;
    clear ax2;
    clear cp_data;
    clear cp_timestamps;
    clear cp_uncertainty;
    clear cruise;
    clear data;
    clear data_to_print;
    clear dataFiles;
    clear eastLon;
    clear end_date;
    clear end_time;
    clear extension;
    clear fid;
    clear fignum;
    clear goodrows;
    clear gps_timestamps;
    clear iData;
    clear j;
    clear lat_data;
    clear lon_data;
    clear northLat;
    clear nRows;
    clear numReplaced;
    clear replaceWith999Index;
    clear sal_data;
    clear sal_timestamps;
    clear sb_fname_ascii;
    clear sb_fname_xls;
    clear sb_hdr;
    clear seabassFileName;
    clear southLat;
    clear start_date;
    clear start_time;
    clear temp_data;
    clear temp_date;
    clear temp_time;
    clear temp_timestamps;
    clear thisData;
    clear thisFileName;
    clear thisPrefix;
    clear thisRow;
    clear timestamps;
    clear timestamps_date_to_print;
    clear timestamps_time_to_print;
    clear title_text;
    clear totalNumReplaced;
    clear wavelengths;
    clear wavelengthToPlot;
    clear westLon;
    clear wlPrefix;
    
    clear a_cut_wl;
    clear a_orig_wl;
    clear ap_std;
    clear ap_std_MatchedUnc;
    clear ap_std_to_print;
    clear c_cut_wl
    clear c_orig_wl;
    clear cp_std;
    clear cp_std_to_print;
    clear cp_std_EndAt750;
    clear data_type;
    clear datamatrix;
    clear date;
    clear fileText;
    clear index;
    clear lat;
    clear lon;
    clear matFileName;
    clear numericCells;
    clear removeIndex;
    clear sal;
    clear std_type;
    clear thisdata;
    clear thisSTD;
    clear time;
    clear timestamp;
    clear titleText;
    clear Wt;
    
    
end;