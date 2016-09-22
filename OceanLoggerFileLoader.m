%OceanLoggerFileLoader - Opens an OceanLogger data file and creates new 
% TemperatureData, SalinityData and FlowData objects
% OceanLoggerFileLoader iterates through a list of files specified in the constructor and uses the
% specified import method to create new TemperatureData, FlowData, and SalinityData objects
%
% Syntax:   obj = OceanLoggerFileLoader(fileNameListIn, importmethodIn, varargin)
%
% Inputs:
%    fileNameListIn - a list of TSG files to import
%    importMethodIn - the import method to use for reading in files
%    varargin       - a set of name/value pairs to specify units for
%    various data objects going to be created
%
%
% Example: 
% olfl = OceanLoggerFileLoader( OceanLoggerFiles, params.TSG_IMPORT_METHOD_NAME, 'temperature', ...
%    params.TEMP_UNITS, 'salinity', params.SAL_UNITS, 'flow', params.FLOW_UNITS );
%
% Other m-files required: TemperatureData, SalinityData, GPSData
% Subfunctions: none
% MAT-files required: none
%
% See also: TemperatureData, SalinityData, FlowData

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% Sept 2016

%------------- BEGIN CODE --------------
classdef OceanLoggerFileLoader
    
    properties
      
        FileNameList     %a list of OceanLogger files to import
        ImportMethodName %the import method to use for reading in files

        
    end
    properties (SetAccess = private, GetAccess = private)
        TemperatureUnits  % Temperature units of measurement
        SalinityUnits     % Salinity units of measurement
        FlowUnits         % Flow units of measurement
    end
    
    methods
        

        %
        % So far this is identical to FlowFileLoader
        %
        function obj = OceanLoggerFileLoader(fileNameListIn, importmethodIn, varargin)
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
                                obj.TemperatureUnits = varargin{i+1};
                            case 'salinity'
                                obj.SalinityUnits = varargin{i+1};
                            case 'flow'
                                obj.FlowUnits = varargin{1+1};
                            otherwise
                                error(['Unexpected option: ' varargin{1}])
                                
                        end
                    end
                end   % end if ~isempty(varargin)
                
            else
                error('Supply an input argument')
            end   % end if nargin > 0
            
        end   % end constructor
             
        

        function dataOut = loadData(obj)   
            %loadData creates a cell array of  data objects for the imported data.
            %It creates a function handle to the import method for this
            %file, allowing the function name to be dynamically set.  The
            %import method currently returns: XXX
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
                [year1,dec_yearday,airtemp1,humidity1,par1,tir1,...
                    airtemp2,humidity2,par2,tir2,baro1,baro2,tstemp,...
                    conductivity,salinity,sound_velocity,chlorophyll,...
                    sampletemp,flowrate,sstemp,trans,sstemp2]...
                    = fh(obj.FileNameList{iFiles});
            
                % Data Processing Section:
                % 
                % Logic is here because this is the object that knows about
                % the specific FILE types -- different files might have
                % different lat/lon type data.  Process here -- GPSData
                % object doesn't need to know format from file
                %

                fullDate = doy2date(dec_yearday, year1);
                
                % ------------------------------------------------------
                % if temperatureData already exists
                if ismember({'td'}, who)
                    
                    % append temperature data
                    td.addData('Temperature', [sstemp,sstemp2], fullDate);
                    
                    % append salinity data
                    sd.addData('Salinity', salinity, fullDate);
                    
                    % append flow Data
                    fd.addData('Flow', flowrate, fullDate);
                else
                    
                    % create new TemperatureData object with the timeseries object
                    td = TemperatureData('Temperature', [sstemp,sstemp2], fullDate, obj.TemperatureUnits);
                                        
                    % create new SalinityData object with the timeseries object
                    sd = SalinityData('Salinity', salinity, fullDate, obj.SalinityUnits); 
                    
                    % create new FlowData object with the timeseries object
                    fd = FlowData('Flow', flowrate, fullDate, obj.FlowUnits);
                    
                end   % end check for existing FlowData object
            end   % end for loop
          
            % create a cell array to hold both
            dataOut = {td, sd, fd};
            
        end   % end loadData function
    end   % end methods block
end   % end classDef
