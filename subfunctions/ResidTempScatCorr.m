function [ap_TSalScatCorr, ap_uncorr_ref, fiterr, deltaT] = ResidTempScatCorr(ap_uncorr, cp, wavelengths, psiT, method)
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
    SLADE = false;
    ROTTGERS = false;
    FLAT = false;
    
    if nargin > 0

        % check input variables             
        if isempty( ap_uncorr ) 
            L.error('ResidTempScatCorr', 'Need ap data');
        end;  
        if isempty( cp )
            L.error('ResidTempScatCorr', 'Need cp data');
        end;
        if isempty( wavelengths )
            L.error('ResidTempScatCorr', 'Need vector of wavelengths');
        end
        if isempty( psiT )
            L.error('ResidTempScatCorr', 'Need psiT');
        end
        if isempty( method )
            L.error('ResidTempScatCorr', 'Need correction method');
        end;
        if strcmpi( method, 'Slade')
            SLADE = true;
        elseif strcmpi( method, 'Rottgers');
            ROTTGERS = true;
        elseif strcmpi( method, 'Flat');
            FLAT = true;
        else
            L.error('ResidTempScatCorr', 'Invalid correction method');
        end;
    
        wavelengths = wavelengths(:)';

        % determine number spectra (number of rows)
        % determine number of wavelengths (number of columns)
        [nSpect, nWavelengths] = size(ap_uncorr);

%         disp('size ap_uncorr')
%         nSpect
%         nWavelengths
        % check data is size we expect;
        if size(ap_uncorr)~=size(cp)
            L.error('ResidTempScatCorr', 'ap,cp size mismatch');
        else
            L.info('ResidTempScatCorr','ap,cp size OK');
        end
        if length(wavelengths)~=nWavelengths
            L.error('ResidTempScatCorr', 'Invalid lambda');
        else
            L.info('ResidTempScatCorr','size wavlengths OK');
        end

        % MOVE THIS TO OUTSIDE FUNCTION, IF WE ARE DOING MORE THAN ONE CORRECTION
        % CHOICE?
%         psiT = getPsiT('Sullivan_etal_2006_instrumentspecific.xls', wavelengths)

        % unchanged from WS/TL code:
        opts = optimset('fminsearch');      
        opts = optimset(opts,'MaxIter',20000000); 
        opts = optimset(opts,'MaxFunEvals',20000); 
        opts = optimset(opts,'TolX',1e-8);
        opts = optimset(opts,'TolFun',1e-8);

        % spectral srange for optimization (710 to 750nm)
        NIR = find( wavelengths >=700 & wavelengths <=750);

        % Reference wl for Slade is 730, Rottgers is 715
        ref730 = find( wavelengths >=730, 1, 'first');
        ref715 = find( wavelengths >=715, 1, 'first');

        ap_TSalScatCorr = nan(size(ap_uncorr));
        deltaT = nan(size(ap_uncorr,1),1);
        ap_uncorr_ref = nan(size(ap_uncorr,1),1);
        fiterr = nan(size(ap_uncorr,1),1);
   
        for k = 1:nSpect
            
            if all(isfinite(ap_uncorr(k,:)))
                
%                 disp('data for k------------------------------------------------')
%                 disp(k)    
                % guess for paramters (beamc at lambda0, beamc slope)
                delT = 0;
                
                % calculate bp for this spectra
                bp = cp(k,:) - ap_uncorr(k,:);
                
                % what is the temperature that makes this look the flattest?
                %
                % this calls f_TS function below
                
                if SLADE
                    
                    ref = find( wavelengths >= 730, 1, 'first');
                        
                    % equation 6 from Slade paper
                    % minimization routine
                    [deltaT(k), fiterr(k)] = fminsearch(@f_TS_Slade, 0, opts, ap_uncorr(k,:), cp(k,:), psiT, NIR, ref);

%                     size(deltaT(k))
%                     size(fiterr(k))
                    
                    % equation 5 from Slade paper
                    % combines temperature and scattering correction    
%                     size(ap_uncorr(k,:))
%                     size(psiT.*deltaT(k))
%                     size( ((ap_uncorr(k,ref)-psiT(ref).*deltaT(k))./bp(ref)).*bp)
                    ap_TSalScatCorr(k,:) = ap_uncorr(k,:) - psiT.*deltaT(k) - ...
                            ((ap_uncorr(k,ref)-psiT(ref).*deltaT(k))./bp(ref)).*bp;

                    % WHAT IS THIS DOING?
                    % return the wavelength that it's using to zero out the
                    % spectra.  Dont' use for Rottgers
                    ap_uncorr_ref(k,1) = ap_TSalScatCorr(k,ref) + (ap_uncorr(k,ref)-psiT(ref).*deltaT(k));
                    
                elseif ROTTGERS
                    
                    ref = find( wavelengths >= 715, 1, 'first');
                    
                    if all(isnan(ref))
                        L.error('ResidTempScatCorr','in Rottgers - ref is nan');
                    end;
                    
                    if all(isnan(NIR))
                        L.error('ResidTempScatCorr','in Rottgers - NIR is nan');
                    else
%                         disp('nir ok')
                    end;
                    if all(isnan(psiT))
                        L.error('ResidTempScatCorr','in Rottgers - psiT is nan');
                    else
%                         disp('psiT ok')
                    end;
                    if all(isnan(cp(k,:)))
                        L.error('ResidTempScatCorr','in Rottgers - cp(k) is nan');
                    else
%                         disp('cp ok')
                    end;
                    if all(isnan(ap_uncorr(k,:)))
                        L.error('ResidTempScatCorr','in Rottgers - ap_uncorr(k,:) is nan');
                    else 
%                         disp('apuncorr ok')
                    end;
                    %EB changed to Slade:
                    %[deltaT(k), fiterr(k)] = fminsearch(@f_TS_Rottgers, 0, opts, ap_uncorr(k,:), cp(k,:), psiT, NIR, ref);    
                    [deltaT(k), fiterr(k)] = fminsearch(@f_TS_Slade, 0, opts, ap_uncorr(k,:), cp(k,:), psiT, NIR, ref);
                    
%                     a715 = 0.212*(ap_uncorr(ref)^1.135);
                    a715 = 0.212*(ap_uncorr(k,ref)^1.135);
                    ec = 0.56;
                    
                    ap_TSalScatCorr(k,:) = ap_uncorr(k,:) - (ap_uncorr(k,ref) - a715).* ...
                        ( (ec^-1*cp(k,:) - ap_uncorr(k,:) ) ./ (ec^-1*cp(k,ref) - a715))...
                        - psiT.*deltaT(k);
                    % test one case
%                     if k == 923
%                         disp('a715')
%                         a715
%                         disp('ec')
%                         ec
%                         disp('ap_TSalScatCorr(k,:)')
%                         ap_TSalScatCorr(k,:)
%                         disp('ap_uncorr(k)')
%                         ap_uncorr(k,:)
%                         disp('(ap_uncorr(k,ref) - a715)')
%                         (ap_uncorr(k,ref) - a715)
%                         disp('(ec^-1*cp(k,:) - ap_uncorr(k,:)')
%                         (ec^-1*cp(k,:) - ap_uncorr(k,:))
%                         disp('(ec^-1*cp(k,ref) - a715)')
%                         (ec^-1*cp(k,ref) - a715)
%                         disp('psiT.*deltaT(k);')
%                         psiT.*deltaT(k)
%                     end;

                elseif FLAT
                    L.debug('ResidTempScatCorr', 'FLAT');

%                     x1 = fminsearch(@f_TS_Flat, [0,0], opts, ap_uncorr(k,:), psiT, NIR);

                        [x1, fiterr(k)] = fminsearch(@f_TS_Flat, [0,0], opts, ap_uncorr(k,:), psiT, NIR);

%                     disp('deltaT')
%                     deltaT(k)
%                     disp('fiterr')
%                     fiterr(k)
                    deltaT(:,:) = x1(1);
                    A(:,:) = x1(2);
                    ap_TSalScatCorr(k,:) = ap_uncorr(k,:) - A - psiT.*deltaT(k);
                    
                else
                        L.error('ResidTempScatCorr', 'Need to run correction');
                end;
                    % deltaT should be same for two corrections...? check?
                    % should be tested
                
%                     % no residual scattering in c-tube, just need deltaT
%                     cp_corr = cp(k,:) - psiT(k,:).*deltaT

                L.debug('ResidTempScatCorr', sprintf('k=%s, deltaT=%s, err=%s', ...
                    num2str(k), num2str(deltaT(k)), num2str(fiterr(k))));

            else
%                 disp('NO data for k------------------------------------------------')
%                 disp(k)  
            end;
        end;   % for loop for k

    else
        L.error('ResidTempScatCorr', 'Need input arguments');
        
    end   % end if nargin > 0
return

function costf = f_TS_Slade(delT, ap, cp, psiT, NIR, ref)

    bp = cp - ap;
    
    costf = sum(abs(   ap(NIR) - psiT(NIR).*delT - ((ap(ref)-psiT(ref).*delT)./bp(ref)).*bp(NIR)    ));

return

% not using right now?
%     function costf = f_TS_Rottgers(delT, ap, cp, psiT, NIR, ref)
% 
%         % 0.212*(ap(ref)^1.135)
% 
%         costf = sum( abs(  ap(NIR) - psiT(NIR).*delT - ...
%             ( (ap(ref) - 0.212*(ap(ref)^1.135)).*...
%             ( ( (cp(NIR)/0.56) - ap(NIR) )./( (cp(ref)/0.56) - 0.212*(ap(ref)^1.135) ) ) )));
% 
%     return
    
function costf = f_TS_Flat( x0, ap, psiT, NIR )
    
    % x0(1) is delT
    % x0(2) is A
    costf = sum( abs( ap(NIR) - x0(2) - psiT(NIR).*x0(1)))
    
return
    