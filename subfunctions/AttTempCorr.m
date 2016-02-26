function [cpCorr, fiterr] = AttTempCorr(cpUncorrIn, wavelengthsIn, psiTIn, methodIn, deltaTIn)
%RESIDTEMPSCATCORR - Performs one of two residual temperature scattering
%corrections:
%   1) Slade et al., (2009)
%   2) Slade and Boss (2015)
%
% Syntax:  [ap_TSalScatCorr, ap_uncorr_ref] = ResidTempScatCorr(ap_uncorr, cp, wavelengths, psiT, method)
% Inputs:
%    ap_uncorr -  a matrix of ap data (uncorrected), rows are spectra,
%                 columns are wavelengths
%    cp - a matrix of cp data, rows are spectra, columns are wavelengths
%    wavelengths - a vector of the wavelengths in ap_uncorr and 
%    psiT - wavelengths from Sullivan spreadsheet
%    method - a string, either "Slade" or "Rottgers"
%
% Outputs:
%    ap_TSalScatCorr - Description
%    ap_uncorr_ref - Description
%
% Example: 
%    [ap_TSalScatCorr_bin_data, ap_uncorr_ref] = ResidTempScatCorr( ...
%       ap_uncorr_bin_data, cp_bin_data, cwl, "Slade");
%    [ap_TSalScatCorr_bin_data, ap_uncorr_ref] = ResidTempScatCorr( ...
%       ap_uncorr_bin_data, cp_bin_data, cwl, "Rottgers");
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Wendy Neary, based on original code by Wayne Slade and/or Tom
% Leeuw
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% Sep 2015; Last revision: 05-09-15

%------------- BEGIN CODE --------------

    L = log4m.getLogger();
    
    % initialize variables;
    WITHA = false;
    WITHOUTA = false;
    
    if nargin > 0

        % check input variables             
        if isempty( cpUncorrIn )
            L.error('AttTempCorr', 'Need uncorrected cp data');
        else
            L.debug('AttTempCorr', sprintf('size of cpUncorrIn: %u %u', size(cpUncorrIn)));
        end;
        if isempty( wavelengthsIn )
            L.error('AttTempCorr', 'Need vector of wavelengths');
        end
        if isempty( psiTIn )
            L.error('AttTempCorr', 'Need psiT');
        end
        if isempty( methodIn )
            L.error('AttTempCorr', 'Need correction method');
        else
            % method not empty -- check:
            if strcmpi( methodIn, 'WITHA')
                WITHA = true;
%                 disp('witha')
                % need deltaT from a correction
                if isempty( deltaTIn )
                    L.error('AttTempCorr', 'Need deltaT passed');
                end;
                
                if all(isnan(deltaTIn))
                    L.error('AttTempCorr', 'deltaT is all nan');
                else
                    L.debug('AttTempCorr', 'deltaT is OK');
                    L.debug('AttTempCorr', sprintf('size of deltaTIn: %u %u', size(deltaTIn)));
                end;
                
            elseif strcmpi( method, 'WITHOUTA');
                WITHOUTA = true;
%                 disp('withouta')
                if isempty( deltaTIn )
                    L.info('AttTempCorr', 'No deltaT passed');
                end;
            else
                L.error('AttTempCorr', 'Invalid correction method');
            end;     
        end;  % end if method is empty
    
        wavelengths = wavelengthsIn(:)';

        % determine number spectra (number of rows)
        % determine number of wavelengths (number of columns)
        [nSpect, nWavelengths] = size(cpUncorrIn);

        % check data is size we expect;
        if length(wavelengths)~=nWavelengths
            L.error('AttTempCorr', 'Invalid lambda');
        end

        % unchanged from WS/TL code:
        % NEED TO REVIEW THESE
        opts = optimset('fminsearch');      
        opts = optimset(opts,'MaxIter',20000000); 
        opts = optimset(opts,'MaxFunEvals',200000); 
        opts = optimset(opts,'TolX',1e-8);
        opts = optimset(opts,'TolFun',1e-8);

        % spectral srange for optimization (710 to 750nm)
        NIR = find( wavelengths >=700 & wavelengths <=750);

        cpCorr = nan(size(cpUncorrIn));
        deltaTOut = nan(size(cpUncorrIn,1),1);
        fiterr = nan(size(cpUncorrIn,1),1);
   
        for k = 1:nSpect
            
%             disp('k------------------------------------------------')
%             disp(k)

            if all(isfinite(cpUncorrIn(k,:)))
             
                if WITHA
             
                    cpCorr(k,:) = cpUncorrIn(k,:) - psiTIn.*deltaTIn(k);
                    
                    
                    L.debug('AttTempCorr', sprintf('WITH A: cpCorr(k)=%s, deltaT=%s, psiT=%s', ...
                        num2str(cpCorr(k)), num2str(deltaTIn(k)), num2str(psiTIn)));
                                       
                elseif WITHOUTA
                    
                    % set up inital guesses

                    % find cp of wavelength at 550 for best guess?
                    delT = 0;
                    gam = 1;
                    
                    % find index of wavelength at 550
                    ref = find( wavelengths >= 550, 1, 'first');
                    amp = cpUncorrIn(ref);
                    x0 =  [amp, gam, delT];
                    
                    
%                     disp('ref')
%                     ref
%                     
%                     disp('cpUncorr(ref) (amp)')
%                     amp
                    
                    if all(isnan(cpUncorrIn(k,:)))
                        L.error('AttTempCorr',sprintf('cpUncorr is nan where k: %u',k));
                    else
                        L.debug('AttTempCorr', 'cpUncorr not nan');
                        
                    end
                    if isnan(psiTIn)
                        L.error('AttTempCorr',sprintf('psiT is nan where k: %u',k));
                    end
                    if isnan(NIR)
                        L.error('AttTempCorr',sprintf('NIR is nan where k: %u',k));
                    end
                    if isnan(wavelengths)
                        L.error('AttTempCorr',sprintf('wavelengths is nan where k: %u',k));
                    end                    
                    
                    if isnan( f_ATT(x0, cpUncorrIn(k,:), psiTIn, NIR, wavelengths))
                        L.error('AttTempCorr','f_ATT returns nan')
                    else
                        L.debug('AttTempCorr','f_ATT doesn''t return nan')
                    end
                        
                    
                    [x1, fiterr(k)] = fminsearch(@f_ATT, x0, opts, cpUncorrIn(k,:), psiTIn, NIR, wavelengths);
                    
                    % x1(1) is amplitude
                    % x1(2) is slope
                    deltaTOut(:,:) = x1(3); % deltaT
              
                    cpCorr(k,:) = cpUncorrIn(k,:) - psiTIn.*deltaTOut(k);
                    
                    L.debug('AttTempCorr', sprintf('k=%s, deltaT=%s, err=%s', ...
                        num2str(k), num2str(deltaTOut(k)), num2str(fiterr(k))));
                    
                else
                        L.error('AttTempCorr', 'Need to run correction');
                end;
            end
        end

    else
        L.error('AttTempCorr', 'Need input arguments');
        
    end   % end if nargin > 0
return

function costf = f_ATT(x0, cp, psiT, NIR, wavelengths)

    amp = x0(1);
    gam = x0(2);
    delT = x0(3);

    costf = sum(abs(   cp(NIR) - amp.*(550./wavelengths(NIR)).^gam - psiT(NIR).*delT ));

return

    