% binSizeIn is a datenum
function  [bin_flags, binned_time, numberBins, binIndexNumbers, ...
                bin_data_mean, bin_sample_size, ...
                bin_data_std] = bin(dataIn, timeIn, binSizeIn, binsToMatchIn )
            
%             obj.L.info('ProcessedData.processFSW','Start of Method');  


            % this is the data we need to bin
            data = dataIn;
            time = timeIn;
%             size(time)
%             disp('time to be binned: end')
%             nandatestr(time(end))
            [~, numberCols] = size(data);
            BIN_SIZE = binSizeIn;
            
          % Binned time stamps
            % convert time datenums to datevec type
            time_to_match_vec = datevec(binsToMatchIn);
            % create a datevec of the first time stamp
            start_time_vec = time_to_match_vec(1,:);
            % set the seconds to zero (just have minute)
            start_time_vec(6) = 0;
            % convert it back to datenum
            start_time_to_match = datenum(start_time_vec);
            
            % create a datevec of the last timestamp
            time_to_match_vec(end,:)
            end_time_vec = time_to_match_vec(end,:);
            % set the seconds to zero (just have minutes)
            end_time_vec(6) = 0;
            % convert it back to datenum
            end_time_to_match = datenum(end_time_vec);
            
%             disp('in bin. end time to match:')
%             datestr(end_time_to_match)
            
            % create binned timestamps, with intervals set by BIN_SIZE
            % parameter
            binned_time = start_time_to_match:BIN_SIZE:end_time_to_match;
            binned_time = binned_time';
            numberBins = numel(binned_time);
            
%             disp('check binned_time created')
%             datestr(binned_time(end))

            % create an index of the bins the data is going IN to
            % i.e. size of original data & time
            binIndexNumbers = zeros(size(time));
            binIndexNumbers(:,:) = NaN;

            % set up flags
            bin_flags = zeros(numberBins, numberCols);
            bin_flags(:,:) = NaN;

            % set up data structures for mean, median, std, 
            bin_data_mean = zeros(numberBins, numberCols);
            bin_data_mean(:,:) = NaN;

            bin_data_std = zeros(numberBins, numberCols);
            bin_data_std(:,:) = NaN;

            bin_sample_size = zeros(numberBins, numberCols);
            bin_sample_size(:,:) = NaN;
            
            % loop through each bin (currently minute)
            for iBin=1:numel(binned_time)

                % get an index of which timestamps belong in this bin
                if iBin < numel(binned_time)
                    % its any bin other than the last one
                    thisBinTimestampIndex = ( time >= binned_time(iBin) ) & ( time < binned_time(iBin+1) );
                else
                    % its the last timestamp
                    thisBinTimestampIndex = ( time >= binned_time(iBin) );
                end

                % Fill in the index number of this bin for future processing
                binIndexNumbers(thisBinTimestampIndex) = iBin;

                % Gather the data for this bin
                thisBinData = data(thisBinTimestampIndex,:);
                [numRows, numCols] = size(thisBinData);

                bin_sample_size(iBin,:) = numRows;

                % Set up the flags for this bin:
                % One flag for this bin for EACH wavelength (1x83)
                % Set all flags to 'good' (1)
                thisBinFlags = bin_flags(iBin,:);
                thisBinFlags(:,:) = 1;  % 1 == good

                % check not ALL of these values are NaN
                if ~all(all(isnan(thisBinData)))
%                     obj.L.debug('ProcessedData.processFSW', sprintf('bin #: %u HAS good data', iBin));
%                     
                    
                    % 3.  MEAN & STD
                    %     use nanmean:  "For matrices X, nanmean(X) is a row vector of
                    %     column means, once NaN values are removed."
                    %     thisBinData is N number of timestamps, with 83 columns.
                    %     We want column means -- one for each column
                    %     Assign to row "i" of bin_data_mean

                    % copy into big matrix for this bin
                    % this is final mean calculation
                    bin_data_mean(iBin,:) = nanmean(thisBinData, 1);
                    bin_data_std(iBin,:) = nanstd(thisBinData, 1);
                    
              
                    % copy in final bin flag?
                    bin_flags(iBin,:) = thisBinFlags;
                else
%                      obj.L.debug('ACData.processFSW', sprintf('Bin %u has no valid data', iBin));
                end   % end check for valid data

            end;    
%         obj.L.info('processFSW','End of Method');           
        end   % end processFSW