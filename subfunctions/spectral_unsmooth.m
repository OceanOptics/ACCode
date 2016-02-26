% THIS NEEDS TO TAKE A OR C

function [abscorr] = spectral_unsmooth( wl, ap, acs)   % ~ is wlunc, apunc

% AC-S "un-smoothing" method, shortened from original spectral_decomp.m
%
% Ron Zaneveld, WET Labs, Inc., 2005
% Ali Chase, University of Maine, 2014
% Wendy Neary, University of Maine, 2015
%
% See the following publication for details on the method:
% Chase, A., et al., Decomposition of in situ particulate absorption
% spectra. Methods in Oceanography (2014), http://dx.doi.org/10.1016/j.mio.2014.02.022
%
% INPUTS:
% wl     -  wavelengths associated with measured particulate absorption
%           spectra
% ap     -  measured particulate absorption spectra (can be either ap or
%           cp)
% acs    -  1 (for AC-S data) or 0 (for a different instrument). This code
%           was designed for use with particulate absorption spectra measured
%           using a WETLabs AC-S instrument, and thus it applies a spectral
%           un-smoothing that is specific to that instrument. To use with
%           other types of absorption data, input "0" and the correction
%           applied for AC-S will be bypassed.
%
% The uncertainty values we use are the standard deviation of one-minute
% binned particulate absorption spectra. If the uncertainty is unknown,
% then all wavelenghts will be treated equally, i.e. no spectral weighting
% will be applied. In this case populate apunc with -999 values.
%
% OUTPUTS:
% abscorr - the unsmoothed particulate spectrum
%
%
% See also: SPECTRAL_DECOMP.M
%
% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% May 2015; Last revision: 16-SEP-15
% Revision Note:
% Line 114 was changed from:
% absspec(maxwavel*10+1:size(wavelength))=0;
% to:
% absspec(maxwavel*10+1:max(size(wavelength)))=absspec(maxwavel*10);%0;


% acs=1;
% wl=[400:2:750];
% ap=exp(-0.014*wl);
% ap=exp(-((wl-550).^2)/100);

%------------- BEGIN CODE --------------

disp('start spectral_unsmooth')
ap_all = ap; %pd.apUncorr;
% nanind = sum(isnan(ap_all),2) > 0;
% ap_all(nanind,:) = [];
% ap_all = ap_all(1,:);

% wl = pd.DeviceFile.aWavelengths;
% acs = 1;

abscorr_out = zeros(size(ap_all));
abscorr_out(:) = NaN;
[numRows,numCols] = size(ap_all);

if acs==1;  % if AC-S

    for iRow=1:numRows   % for each timestamp
        ap = ap_all(iRow,:);
        if all(isnan(ap))
            % do nothing, all nan data, just copy 
            abscorr_out(iRow,:) = ap_all(iRow,:);
%             disp('all nan')
        else
    %------------- BEGIN LOOP THROUGH TIMESTAMPS ----------------
    % disp('should be data')
            apsize=size(ap);
            if apsize(1)<apsize(2);
                ap=ap';
            end

            wlsize=size(wl);
            if wlsize(1)<wlsize(2);
                wl=wl';
            end

            % Set up filter factors at every 0.1 nm from 1 to 799 nm, with center
            % wavelength at centwavel (i.e. at the data wavelengths)
            wavelength=350:.1:799; % Thus index of 1 nm = 10; 356 nm= 3560;
            clear filtfunc
            SIG1=(-9.845*10^-8.*wl.^3  +1.639*10^-4*wl.^2- 7.849*10^-2*wl + 25.24)/2.3547 ;
            for i=1:max(size(wl));
                for jkl=1:max(size(wavelength));
                    % First term normalizes area under the curve to 1:
                    filtfunc(jkl,i)=(1/(sqrt(2*pi)*SIG1(i)))*exp(-0.5*((wavelength(jkl)-wl(i))/SIG1(i)).^2); 
                end
            end

            error=10^-2;    % stop iterating when error in calculated ap is less than this from measured ap.
            numbits=10 ;    % number of iterations

            % Convolve the measurement with the fiter factors add the difference to
            % the measured spectrum to get the first corrected spectrum and iterate
            % to find the spectrum that when convolved with the filters, will give the
            % measured spectrum. This is the corrected absorption spectrum "abscorr".

            minwavel=min(wl);
            maxwavel=max(wl);

            abscorr = zeros(size(ap));
            abscorr(:) = NaN;

            for ii=numbits;
                xixi=minwavel:.1:maxwavel;% The range of centwavel is 0.1 nm.
                yiyi=spline(wl,ap,xixi); % Spline the measured data to every 0.1 nm.
                % We need data from 350nm to 799 nm to multiply by filtfac.
                absspec=zeros(size(wavelength));
                absspec((minwavel-350)*10+1:(maxwavel-350)*10+1)=yiyi;
                absspec(1:(minwavel-350)*10)=interp1(wl,ap,wavelength(1:(minwavel-350)*10),'linear','extrap');
                absspec((maxwavel-350)*10+2:length(wavelength))=absspec((maxwavel-350)*10+1);%0;
                aspecprime=absspec';

                clear meassignal6
                
                meassignal6=0.1*filtfunc'*aspecprime;
                
                abscorr=ap-meassignal6+ap;
                if max((ap-meassignal6)/ap)<=error
                    break
                end     

            end  % for ii=numbits

            abscorr_out(iRow,:) = abscorr(:,:)';

        end  % end check for nans
        disp('-');
            %--------------------------end code

    end  % end loop through timestamps 
    %---------------------------------- end loop through timestamps
else  % not an ac-s

    abscorr_out=ap;
%     abscorr_out(iRow,:) = ap_all(iRow,:)';
end

% set return variable
abscorr = abscorr_out;

%end   % end function

%------------- END OF CODE --------------



