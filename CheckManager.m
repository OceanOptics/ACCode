%% CheckManager - SCRIPT which inputs SeaBASS data files and plots them
% for visual inspection
%
%
% Requires: the .mat file saved in ProcessingManager, if not in memory
%           the .mat file saved in PreProcessingManager, if not in memory
%           SeaBASS text data files
%
% Outputs:  Plots for ap_rottgers, ap_slade, cp
%
% Other m-files required: importSeaBASS2 
% Subfunctions: checkTimestamps
% MAT-files required: "acsPROC_(YEAR)_(DAY)" saved by ProcessingManager,
%  "acsPREPROC_(YEAR)_(DAY)" saved by PreProcessingManager
%
% See also: ProcessedData, importSeaBASS2
% 
% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% Jan 2016; Last revision: 25-Apr-16
% 25-Apr-16: Added new file format that includes std and bin_count in
% SeaBASS file formats
% 
%----------------------------- BEGIN CODE ---------------------------------
%%Load data file from disk

if params.RUN.LOAD_CHECK_DATA_FROM_DISK

    matFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
        strcat('acsPROC', '_', num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY)));
    load(matFileName);

    matFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
        strcat('acsPREPROC', '_', num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY)));
    load(matFileName);  
   
end;

%get wavelengths
wavelengths = pd.var.ap.L8.wavelengths_slade;
sb_fname_ascii = [params.OUTPUT.SEABASS_FILE_PREFIX ...
    num2str(params.INGEST.YEAR) '_' num2str(params.INGEST.YEAR_DAY)];

dataFiles = {'ap', 'cp'}; 
fileText = {'ap - corrected using Slade', 'cp'};
wlPrefix = {'ap', 'cp'};               

%%
for iData = 1:length(dataFiles)

    thisFileName = dataFiles{iData}; %;
    thisPrefix = wlPrefix{iData}; %;
    
    extension = strcat(thisFileName, '.txt');
    seabassFileName = fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
        strcat(sb_fname_ascii, extension));

    datamatrix = importSeaBASS3(seabassFileName);

    date = datamatrix(:,1);
    time = datamatrix(:,2);
    lat = datamatrix(:,3);
    lon = datamatrix(:,4);
    Wt = datamatrix(:,5);
    sal = datamatrix(:,6);
%     numericCells = datamatrix(:,7:end-1);

    numWavelengths = length(wavelengths);  
    endColumn = 7 + numWavelengths - 1;
    numericCells = datamatrix(:,7:endColumn);
    
%     thisdata = cell2mat(numericCells);
% changed 4/25/16
    thisdata=cellfun(@str2num,numericCells);
    index = (thisdata(:,:) == -9999);
    thisdata(index) = NaN;

    timestamp = datenum(strcat(date,time),'yyyymmddHH:MM:SS');
    fignum = iData+50;
    
    titleText = fileText{iData};

    figure(fignum)
    hold on;
    grid on;
    plot(wavelengths, thisdata)
    xlabel('Wavelength')
    
    ylabel(thisFileName)
    title_text = sprintf('SeaBASS data\n %s vs. Wavelength\n Cruise: %s Leg: %s \n%u-%s-%u (Yearday: %u)', ...
        titleText, params.INGEST.CRUISE, params.INGEST.CRUISE_LEG, params.INGEST.DAY, ...
        params.INGEST.MONTH_TEXT, params.INGEST.YEAR, params.INGEST.YEAR_DAY);
    title(title_text, 'fontsize', 12);
    
    saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
        strcat(num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY), ...
        '_CHECK_', thisFileName)));
    saveas(gcf,  fullfile(params.INGEST.DATA_OUTPUT_DIRECTORY, ...
        strcat(num2str(params.INGEST.YEAR), '_', num2str(params.INGEST.YEAR_DAY), ...
        '_CHECK_', thisFileName, '.jpg')));  

    
end;

if params.INGEST.CLEAR_VARS
%     clear pd;
end;