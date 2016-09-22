%SEATEXFileLoader - Opens a Seatex data file and creates new GPSData object
% SEATEXFileLoader iterates through a list of files specified in the constructor and uses the
% specified import method to create new GPSData object
%
% Syntax:   obj = SEATEXFileLoader(fileNameListIn, importmethodIn, varargin)
%
% Inputs:
%    fileNameListIn - a list of TSG files to import
%    importMethodIn - the import method to use for reading in files
%    varargin       - a set of name/value pairs to specify units for
%    various data objects going to be created
%
%
% Example: 
% tfl = SEATEXFileLoader( SEATEXFiles, params.SEATEX_IMPORT_METHOD_NAME,...
%    'gps', params.GPS_UNITS );
%
% Other m-files required: TemperatureData, SalinityData, GPSData
% Subfunctions: none
% MAT-files required: none
%
% See also: GPSData

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% Sep 2106

%------------- BEGIN CODE --------------
classdef SEATEXFileLoader
    
    properties
      
        FileNameList     %a list of SeaTex files to import
        ImportMethodName %the import method to use for reading in files

        
    end
    properties (SetAccess = private, GetAccess = private)
        GPSUnits         % GPS units of measurement
    end
    
    methods
        

        %
        % So far this is identical to FlowFileLoader
        %
        function obj = SEATEXFileLoader(fileNameListIn, importmethodIn, varargin)
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
                            case 'gps'
                                obj.GPSUnits = varargin{1+1};
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
            %import method currently returns: XXXX
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
                [year1, dec_yearday, lat, lon, gq, svc, hdop, altitude] = ...
                    fh(obj.FileNameList{iFiles});
            
                % Data Processing Section:
                % 
                % Logic is here because this is the object that knows about
                % the specific FILE types -- different files might have
                % different lat/lon type data.  Process here -- GPSData
                % object doesn't need to know format from file
                %
                % date format here hardcoded -- change?
                fullDate = doy2date(dec_yearday, year1);
                
                % -------------------------------------------------------
                % Put latitude and longitude into one matrix
                LatLon = [lat,lon];
                
                % ------------------------------------------------------
                % if FlowData already exists
                if ismember({'gd'}, who)
                  
                    % append GPS Data
                    gd.addData('GPS', LatLon, fullDate);
                else
                    
                    % create new GPSData object with the timeseries object
                    gd = GPSData('LatAndLon', LatLon, fullDate, obj.GPSUnits);
                    
                end   % end check for existing FlowData object
            end   % end for loop
          
            % create a cell array to hold both
            dataOut = {gd};
            
        end   % end loadData function
    end   % end methods block
end   % end classDef
