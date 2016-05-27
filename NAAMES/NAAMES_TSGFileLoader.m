%TSGFileLoader - Opens a TSG data file and creates new TemperatureData, GPSData, and SalinityData objects
% TSGFileLoader iterates through a list of files specified in the constructor and uses the
% specified import method to create new TemperatureData, GPSData, and SalinityData objects
%
% Syntax:   obj = TSGFileLoader(fileNameListIn, importmethodIn, varargin)
%
% Inputs:
%    fileNameListIn - a list of TSG files to import
%    importMethodIn - the import method to use for reading in files
%    varargin       - a set of name/value pairs to specify units for
%    various data objects going to be created
%
%
% Example: 
% tfl = TSGFileLoader( tsgFiles, params.TSG_IMPORT_METHOD_NAME, 'temperature', ...
%    params.TEMP_UNITS, 'salinity', params.SAL_UNITS, 'gps', params.GPS_UNITS );
%
% Other m-files required: TemperatureData, SalinityData, GPSData
% Subfunctions: none
% MAT-files required: none
%
% See also: TemperatureData, SalinityData, GPSData

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% May 2015; Last revision: 13-08-15

%------------- BEGIN CODE --------------
classdef NAAMES_TSGFileLoader
    
    properties
      
        FileNameList     %a list of TSG files to import
        ImportMethodName %the import method to use for reading in files

        
    end
    properties (SetAccess = private, GetAccess = private)
        TemperatureUnits % Temperature units of measurement
        SalinityUnits    % Salinity units of measurement
        GPSUnits         % GPS units of measurement
        
%         % ADDED NEW LINES HERE:
%         FLRUnits
    end
    
    methods
        

        %
        % So far this is identical to FlowFileLoader
        %
        function obj = NAAMES_TSGFileLoader(fileNameListIn, importmethodIn, varargin)
            % constructor -- takes a list of files and sets object properties 

            % check incoming arguments
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
                
              
                if ~isempty(varargin)
                    for i = 1:2:length(varargin)
                        val = lower(varargin{i});
                        switch val

                            case 'temperature'
                                disp('flow')
                                obj.TemperatureUnits = varargin{i+1};
                            case 'salinity'
                                disp('valve')
                                obj.SalinityUnits = varargin{i+1};
                            case 'gps'
                                obj.GPSUnits = varargin{1+1};
                                %ADDED NEW LINES HERE:
%                             case 'flr'
%                                 obj.FLRUnits = varargin{1+1};
                            otherwise
                                error(['Unexpected option: ' varargin{1}])
                                
                        end
                    end
                end   % end if ~isempty(varargin)
                
            else
                error('Supply an input argument')
            end   % end if nargin > 0
            
        end   % end constructor
             
        

        function dataOut = loadData(obj)   %tscOut = loadData(obj)
            %loadData creates a cell array of  data objects for the imported data.
            %It creates a function handle to the import method for this
            %file, allowing the function name to be dynamically set.  The
            %import method currently returns: [TemperatureBoat, TemperatureInstr, Salinity, LatDegrees, ...
            %        LatMinutes, LatSign, LonDegrees, LonMinutes, LonSign, ...
            %        Time1, Date1]
            %
            % SYNOPSIS dataOut = loadData(obj)
            % INPUT obj      - this object
            % OUTPUT dataOut - a cell array of TemperatureData,
            % SalinityData and GPSDataobjects
            % 
                   
           
            nFiles = length(obj.FileNameList);
            for iFiles = 1:nFiles;
                
                
                % open the data file and retrieve data from it
                fh = str2func(obj.ImportMethodName);
                [DATE_GMT,TIME_GMT,Dec_LAT,Dec_LON,SBE45S,SBE48T,FLR] ...
                    = fh(obj.FileNameList{iFiles});
            
                % Data Processing Section:
                % 
                % Logic is here because this is the object that knows about
                % the specific FILE types -- different files might have
                % different lat/lon type data.  Process here -- GPSData
                % object doesn't need to know format from file
                %
                % date format here hardcoded -- change?

                fullDate = datenum(strcat(DATE_GMT,TIME_GMT),  'yyyy/mm/ddHH:MM:SS.FFF');
                
                % -------------------------------------------------------
                % Rename Lat & Lon:
                Latitude = Dec_LAT;
                Longitude = Dec_LON;
                
                % Not doing this for NAAMES:
%                 % Set sign appropriate for N/S and E/W
%                 Latitude(LatSign=='S') = -Latitude(LatSign=='S');
%                 Longitude(LonSign=='W') = -Latitude(LonSign=='W');
                
                % Put latitude and longitude into one matrix
                LatLon = [Latitude, Longitude];
                
                % Rename Salinity:
                Salinity = SBE45S;
                TemperatureBoat = SBE48T;
                TemperatureInstr = SBE48T;
                
                %ADDED NEW LINE HERE:
                % Rename FLR
%                 FLR = FLR;
                % ------------------------------------------------------
                % if TSGData already exists
                if ismember({'td'}, who)
                    
                    %append data
                    % Use TemperatureInstr as the temperature measurement
                    td.addData('Temperature', [TemperatureInstr, TemperatureBoat], fullDate);
                    sd.addData('Salinity', Salinity, fullDate);
                    
                    % append GPS Data
                    gd.addData('GPS', LatLon, fullDate);
                    
%                     % ADDED NEW LINE HERE:
%                     % append FLR Data
%                     fd.addData('FLR', FLR, fullDate);
                else
                    
                    % create new TemperatureData object with the timeseries object
                    td = TemperatureData('Temperature', [TemperatureInstr, TemperatureBoat], fullDate, obj.TemperatureUnits);
                                        
                    % create new SalinityData object with the timeseries object
                    sd = SalinityData('Salinity', Salinity, fullDate, obj.SalinityUnits); 
                    
                    % create new GPSData object with the timeseries object
                    gd = GPSData('LatAndLon', LatLon, fullDate, obj.GPSUnits);
                    
%                     % ADDED NEW LINE HERE:
%                     % create new FLRData object with the timeseries object
%                     fd = FLRData('FLR', FLR, fullDate, obj.FLRUnits);
                    
                end   % end check for existing FlowData object
            end   % end for loop
          
            % create a cell array to hold both
            % ADDED 'fd' HERE:
%             dataOut = {td, sd, gd, fd};
            dataOut = {td, sd, gd};
            
        end   % end loadData function
    end   % end methods block
end   % end classDef
