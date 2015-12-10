%% testFlowImport
%
% This script tests the object FlowData
% --------------------------------------------------------------------------
%%
% get rid of all old objects.  use delete to get rid of handle class
clear all;

params = readIngestParameters();
% generate list of flow files
   
flowcurr = sprintf(params.FLOW_FILE_FORMAT, params.YEAR, params.YEAR_DAY);
flownext = sprintf(params.FLOW_FILE_FORMAT, params.YEAR, params.YEAR_DAY+1);
flowprev = sprintf(params.FLOW_FILE_FORMAT, params.YEAR, params.YEAR_DAY-1);
flowFiles = getFilesNextPrev(params.FLOW_DIRECTORY, flowcurr, flownext, flowprev);

ffl = FlowFileLoader(flowFiles, params.FLOW_IMPORT_METHOD_NAME, 'flow', params.FLOW_UNITS, 'valve', params.VALVE_UNITS );
flowData = ffl.loadData();

flow = flowData{1};
%%
% ingest variables and methods
% display variables of new object
flow.Name
flow.Type
flow.DataObject
disp('check quality')
flow.DataObject.QualityInfo

% test methods are being called properly
flow.getInfo()
figure(1)
scatter(flow.DataObject.Time, flow.DataObject.Data, 'k')
dynamicDateTicks();
title('Flow')

valve = flowData{2};
figure(2)
scatter(valve.DataObject.Time, valve.DataObject.Data)
dynamicDateTicks();
title('Valve State')







