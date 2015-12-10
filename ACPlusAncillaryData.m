% ACPLUSANCILLARYDATA is a class to contain the data for AC Data
% Processing.  It can contain a variable number of data types, but must
% contain the following:
% PROPERTIES
%   aData - absorption data as measured by the ac meter
%   cData - attenuation data as measured by the ac meter
%   DeviceFile - information about the AC meter actually deployed on the
%   cruise this ac data comes from
%   GPSData - latitude and longitude data
%   TemperatureData - temperature measurements (could be more than one)
%   SalinityData - salinity measurements
% VARIABLE PROPERTIES
%   FlowData - data taken from a flow metere that determines the amount of
%   flow through the system at a given point
%   ValveData - data about whether the filtered/unfiltered valve is
%   open/closed
%   XYZData - additional data nestuffeded for other data processing

classdef ACPlusAncillaryData
    
    properties (Access = private)
        mandatoryargin = 6;
    end
    properties
        % mandatory properties:
        DeviceFile
        IngestParameters
        
        aData
        cData
        GPSData
        TemperatureData
        SalinityData
        
        % optional properties:
        FlowData
        ValveData
    end
    methods
        function obj = ACPlusAncillaryData( DeviceFileIn, DataTypesIn, IngestParamsIn )
            nargin
            if nargin > 0 
                % set device file
                if isa( DeviceFileIn, 'ACDeviceFile' )
                    obj.DeviceFile = DeviceFileIn;
                else
                    error('Supply DeviceFile object')
                end
                
                if isa( IngestParamsIn, 'struct' )
                    obj.IngestParameters = IngestParamsIn;
                else
                    error('Supply Ingest Params')
                end
                
                % go through other args
                nTypes = length(DataTypesIn);
                for k = 1:nTypes
                    
                    dataName = class( DataTypesIn{k} );
                    
                    switch dataName
                        case 'ACData'
                            if strtok( DataTypesIn{k}.Name, ':' ) == 'a'
                                obj.aData = DataTypesIn{k};
                            elseif strtok( DataTypesIn{k}.Name, ':' ) == 'c'
                                obj.cData = DataTypesIn{k};
                            end
                        case 'GPSData'
                            obj.GPSData = DataTypesIn{k};
                        case 'TemperatureData'
                            obj.TemperatureData = DataTypesIn{k};
                        case 'SalinityData'
                            obj.SalinityData = DataTypesIn{k};
                        case 'FlowData'
                            obj.FlowData = DataTypesIn{k};
                        case 'ValveData'
                            obj.ValveData = DataTypesIn{k};
                    end   % end switch statement
                    
                end   % end loop through varargin
                
            end   % if nargin > 2
            
        end % end function def
        
        
    end   % end methods block

end   % end classdef block
            
            
            