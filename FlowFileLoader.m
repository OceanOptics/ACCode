%FlowFileLoader - Open a Flow data file and create a new ValveData & FlowData objects
%FlowFileLoader iterates through a list of files specified in the 
%constructor and uses the specified import method to create a ValveData and FlowData object.
%Import method must return [Date1, Time1, ValveState, Flow1]
%
% Syntax:  obj = FlowFileLoader(fileNameListIn, importmethodIn, varargin)
%
% Inputs:
%    fileNameListIn - a list of full file names of the flow files
%    importMethodIn - the name of the import script to use to read files
%    Name           - the name to specify which units are being given
%    Value          - the units for a given name
%
% Example: 
%    ffl = FlowFileLoader(flowFiles, params.FLOW_IMPORT_METHOD_NAME, 'flow', params.FLOW_UNITS, 'valve', params.VALVE_UNITS );
%
% Other m-files required: FlowData, log4m
% Subfunctions: none
% MAT-files required: none
%
% See also: FlowData,  IngestManager

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% May 2015; Last revision: 13-08-15

%------------- BEGIN CODE --------------
classdef FlowFileLoader
    
    properties
        FileNameList       % list of files to read in
        ImportMethodName   % name of import method to use to read files
        
    end
    properties (SetAccess = private, GetAccess = private)
        FlowUnits
        ValveUnits
        L   % Logger
    end
    methods
        

        function obj = FlowFileLoader(fileNameListIn, importmethodIn, varargin)
        %FlowFileLoader creates the FlowFileLoader object.
        %
        % SYNOPSIS obj = FlowFileLoader(fileNameListIn, importmethodIn, varargin)
        % INPUT  
        %    fileNameListIn - a list of full file names of the flow files
        %    importMethodIn - the name of the import script to use to read files
        %    Name           - the name to specify which units are being given
        %    Value          - the units for a given name
        % OUTPUT obj -  this object
        % 
         
            
            obj.L = log4m.getLogger();
            obj.L.debug('FlowFileLoader.FlowFileLoader()','Created object');
            
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
                
                % parse varargin to get units
                if ~isempty(varargin)
                    for i = 1:2:length(varargin)
                        val = lower(varargin{i});
                        switch val

                            case 'flow'
                                obj.FlowUnits = varargin{i+1};
                            case 'valve'
                                obj.ValveUnits = varargin{i+1};
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
            %import method must return: [Date1, Time1, ValveState, Flow1]
            %
            % SYNOPSIS dataOut = loadData(obj)
            % INPUT obj      - this object
            % OUTPUT dataOut - a cell array of FlowData and ValveData objects
            % 
                   
            % iterate through list of files
            % if no FlowData object, create it.
            % if FlowData object, append to it.
            
            nFiles = length(obj.FileNameList);
            for iFiles = 1:nFiles;
                
                obj.L.debug('FlowFileLoader.loadData', sprintf('current iteration: %u', iFiles));
                
                % open the data file and retrieve data from it
                fh = str2func(obj.ImportMethodName);
                [Date1, Time1, ValveState, Flow] = fh(obj.FileNameList{iFiles});
            
%                 Date1(1)
%                 Time1(1)
                
                % date format here hardcoded -- change?
                xx = strcat(Date1, Time1);
                fullDate = datenum(xx, 'yyyy-mm-ddHH:MM:SS');
                
                % if FlowData already exists
                % call this 'Flow1' and 'Valve1' because this is the first
                % flow file used with this data. If we load flow data from
                % another flow file, that could be 'Flow2'
                if ismember({'fd'}, who)
                    
                    %append data
                    fd.addData('Flow1', Flow, fullDate);
                    vd.addData('Valve1', ValveState, fullDate);
                   
                else
%                     Flow
%                     fullDate
                    % create new FlowData object with the timeseries object
                    fd = FlowData('Flow1', Flow, fullDate, obj.FlowUnits);
                                        
                    % create new ValveData object with the timeseries object
                    vd = ValveData('Valve1', ValveState, fullDate, obj.ValveUnits); 
                    
                end   % end check for existing FlowData object
            end   % end for loop
         
            % create a cell array to hold both
            dataOut = {fd, vd};
            
        end   % end loadData function
    end   % end methods block
end   % end classDef
