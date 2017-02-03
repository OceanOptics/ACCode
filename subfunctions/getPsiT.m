function [ psiT ] = getPsiT( fileNameIn, wavelengthsIn )
%FUNCTION_NAME - One line description of what the function or script performs (H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [ psiT ] = getPsiT( fileNameIn, wavelengthsIn )
%
% Inputs:
%    fileNameIn - An excel spreadsheet of temperature corrections  -
%    currently ("Sullivan_etal_2006_instrumentspecific.xls)
%    wavelengthsIn - A list of wavelengths to look up corrections for
%
% Outputs:
%    psiT - a temperature correction for each wavelength
%
% Example: 
%    psiT = getPsiT('Sullivan_etal_2006_instrumentspecific.xls', cwl)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: RESIDTEMPSCATCORR,  

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% Sep 2015; Last revision: 05-SEP-15

%------------- BEGIN CODE --------------
    L = log4m.getLogger();
    if nargin > 0

        % check there is a file list              
        if ~isempty( fileNameIn ) 
            fileName = fileNameIn;
        else
            L.error('getPsiT', 'Need file name');
        end  
        
        % check there are wavelengths
        if ~isempty( wavelengthsIn )
            wavelengths = wavelengthsIn;
        else
            L.error('getPsiT', 'Need vector of wavelengths');
        end

        % open file
        tmp = xlsread('Sullivan_etal_2006_instrumentspecific.xls');
        if ~isempty(tmp)
            % added extrap 11/4/16
            psiT = interp1(tmp(:,1),tmp(:,2), wavelengths, 'linear', 'extrap');
        else
            L.error('getPsiT', 'Xlsread returned null');
        end;
        
    else   % no args
        L.error('getPsiT', 'Need input arguments');
    end;

end

