function [ goodStartsOut, goodEndsOut, startFlagsOut, endFlagsOut ] = ...
    checkTransitionTimes( rawStartsIn, rawEndsIn, startFlagsIn, endFlagsIn, timestampsIn ) 
%CHECKTRANSTIONTIMES - This function matches a set of start and end
%transition times.
% If data starts with a transition end, it adds the first timestamp from
% timestampsIn as that transition start.  Similarly, if data ends with a 
% transition start, it adds the last timestamp as the end for that start.
% If a transition start or end is missing, NaN is entered in that list and
% the flag is set to '3'.
% The output flags use 1 == good, 3 == suspect.  

% Any original flags are
% copied over, but if a transition is bad, that flag will be set to
% suspect.
% 
%
% Syntax:  [ goodStartsOut, goodEndsOut, startFlagsOut, endFlagsOut ] = ...
%    checkTransitionTimes( rawStartsIn, rawEndsIn, startFlagsIn, endFlagsIn, timestampsIn )
%
% Inputs:
%    rawStartsIn  - a list of datenums of identified transition starts
%    rawEndsIn    - a list of datenums of identified transition ends
%    startFlagsIn - the original flags associated with these start times
%    endFlagsIn   - the original flags associated with these end times
%    timestampsIn - the original timestamps these transtion times came from
%
% Outputs:
%    goodStartsOut - a list of transition start times, matched to 
%           goodEndsOut.  Contains NaNs.
%    goodEndsOut   - a list of transitions end times, matched to
%           goodStartsOut.  Contains NaNs.
%    startFlagsOut - the flags associated with goodStartsOut
%    endFlagsOut   - the flags associated with goodEndsOut
%
% Example: 
%    [ goodStarts, goodEnds, startFlags, endFlags ] = ...
%    checkTransitionTimes( starts, ends, startflags, endflags, timestamps ); 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Wendy Neary
% MISCLab, University of Maine
% email address: wendy.neary@maine.edu 
% Website: http://misclab.umeoce.maine.edu/index.php
% Nov 2015; Last revision: 05-11-15

%------------- BEGIN CODE --------------


    L = log4m.getLogger();

    % assign input variables
    rawStartTimes = rawStartsIn;
    rawEndTimes = rawEndsIn;
    timestamps = timestampsIn;

    % these are incoming flags that might have been set for a different reason
    startFlags = startFlagsIn;
    endFlags = endFlagsIn;

    startsLength = length(rawStartTimes);
    endsLength = length(rawEndTimes);

    % which is longer?
    if startsLength >= endsLength
        L.debug('matchTransitionStartsEnds', ...
            sprintf('more or equal starts (%u) than ends: %u', startsLength, endsLength));
        maxLength = startsLength;
    else
        L.debug('matchTransitionStartsEnds', ...
            sprintf('fewer starts (%u) than ends: %u', startsLength, endsLength));
        maxLength = endsLength + 1;
    end;

    % ------------------------------------------------------------
    % 1.  Loop through and match starts with ends
    %---------------------------------------------------------------

    % calculate number of ends before first beginning
    % could be more than one if false ends
    tempidx = rawEndTimes < rawStartTimes(1);
    numEndsWithNoStarts = sum(tempidx);
    L.debug('matchTransitionStartsEnds', ...      
       sprintf('numEndsWithNoStarts: %u', numEndsWithNoStarts));

    % we have one end with no start at the beginning
    if numEndsWithNoStarts == 1

       % if we have ONE end with no starts, increase the array of
       % starts by ONE
       tempStarts = zeros( length(rawStartTimes) + 1, 1);
       tempStarts(:,:) = NaN;

       % we have some end times before our first start time
       tempStarts(2:end) = rawStartTimes(:,:);

       % assign tempStarts(1) first timestamp in data file
       tempStarts(1) = timestamps(1);

       % do same for flags
       tempStartFlags = zeros(numEndsWithNoStarts+length(rawStartTimes), 1);
       tempStartFlags(:,:) = NaN;

       L.debug('checkTransitionTimes', ...
           sprintf('size numEnds+1: %u', length(tempStartFlags(numEndsWithNoStarts+1:end))));

       % copy in any flags already set by previous code
       tempStartFlags(numEndsWithNoStarts+1:end) = startFlags;

       L.debug('checkTransitionTimes', ...
           sprintf('number good flags: %u', sum(tempStartFlags==1)));

       % mark the nan flags as good
       % could be marked as something else, but with a note that it was
       % changed?
       tempStartFlags(1:numEndsWithNoStarts) = 1;
       L.debug('checkTransitionTimes', ...
           sprintf('number good flags after fixing nans: %u', sum(tempStartFlags==1)));

       % set up goodStarts and goodEnds to copy data into
       goodStarts = zeros(numEndsWithNoStarts+length(rawStartTimes), 1);
       goodStarts(:,:) = NaN;
       % create same size flags
       goodStartFlags = goodStarts;

       goodEnds = zeros(numEndsWithNoStarts+length(rawStartTimes), 1);
       goodEnds(:,:) = NaN; 
       % create same size flags
       goodEndFlags = goodEnds;
       L.debug('checkTransitionTimes', ...
           sprintf('size goodEndFlags: %u', size(goodEndFlags)));
       % set up end Flags: same length as ends, copy over
       % changed 12/1/15
       tempEndFlags = zeros(length(goodEndFlags), 1);
              L.debug('checkTransitionTimes', ...
           sprintf('size tempEndFlags: %u', size(tempEndFlags)));
       tempEndFlags(:,:) = NaN;
       tempEndFlags(1:length(endFlags),:) = endFlags(:,:);
              L.debug('checkTransitionTimes', ...
           sprintf('size tempEndFlags: %u', size(tempEndFlags)));
    % no end at start without beginning   
    elseif numEndsWithNoStarts == 0

       % no ends before first start time
       tempStarts = zeros(maxLength, 1);
       tempStarts(:,:) = NaN;
       tempStarts = rawStartTimes;

       tempStartFlags = startFlags; 
       

       goodStarts = zeros(maxLength, 1);
       goodStarts(:,:) = NaN;
       goodStartFlags = goodStarts;

       goodEnds = zeros(maxLength, 1);
       goodEnds(:,:) = NaN; 
       goodEndFlags = goodEnds;
%        disp('***********************************************')
       tempEndFlags = zeros(length(goodEndFlags),1);
       tempEndFlags(:,:) = NaN;
       tempEndFlags = endFlags;

    else   % numEndsWithNoStarts > 1

       L.error('checkTransitionTimes', ...
           'there should not be more than one end without a start');
    end;

    % for each good flagged start time, find the corresponding end
    % time and copy it into a new array (don't worry about it
    % changing size on each iteration)
    % look for end marker at end of period
    L.debug('checkTransitionTimes', ...
       sprintf('# tempStarts %u', length(tempStarts)));
    L.debug('checkTransitionTimes', ...
        sprintf('# good flags %u', sum(tempStartFlags==1)));
    L.debug('checkTransitionTimes', ...
        sprintf('size goodStarts %u', size(goodStarts)));
    L.debug('checkTransitionTimes', ...
        sprintf('size goodEnds %u', size(goodEnds)));    

    numTempStarts = length(tempStarts);

    % need double loop here - one to loop through the temp starts/ends;
    % two to loop through good starts ends;
    numGoodStarts = length(goodStarts);
    numGoodEnds = length(goodEnds);
    iGoodStart = 1;
    iGoodEnd = 1;
    
    for iTempStart = 1:numTempStarts;
    L.debug('checkTransitionTimes', ...
        sprintf('looping through tempStarts: %u:', iTempStart));
                
        % if this start is not nan - why would a start be nan?
        if ~isnan(tempStarts(iTempStart))
                   
           L.debug('checkTransitionTimes', ...
              sprintf('current one not nan:-------------- %s', datestr(tempStarts(iTempStart))));

          % find the number of possible ends
            if iTempStart < numTempStarts
                % Look for times greater than current time, but
                % less than next time
                possibleEndsIndex = rawEndTimes > tempStarts(iTempStart) &  rawEndTimes < tempStarts(iTempStart + 1);
                numEndsFound = sum(possibleEndsIndex);
            else
                % last one
                % look for times greater than this current start 
                possibleEndsIndex = rawEndTimes > tempStarts(iTempStart);
                numEndsFound = sum(possibleEndsIndex);
            end;
                        
                    
           if numEndsFound == 1
               % we have one exact match
              L.debug('checkTransitionTimes', ...
                  'we have exactly one match');
              
              % <------- COPY IN DATA 1/3 -------------------------->
               goodStarts(iGoodStart) = tempStarts(iTempStart);
               goodStartFlags(iGoodStart) = tempStartFlags(iTempStart);

               goodEnds(iGoodEnd) = rawEndTimes(possibleEndsIndex);
               % changed 11/17 from:
%                goodEndFlags(iGoodEnd) = 1; % interpolated?
                % to: (copy over OLD flag)
%                 size(goodEndFlags)
%                 size(tempEndFlags)
%                 iGoodEnd;
%               iTempStart
              if ~isnan(tempEndFlags(iTempStart))
                  
                 goodEndFlags(iGoodEnd) = tempEndFlags(iTempStart); 
              else
                   L.debug('checkTransitionTimes', ...    
                       sprintf('last tempEndFlag is nan, copying in 1'));
                   goodEndFlags(iGoodEnd) = 1;
              end;
                 
               L.debug('checkTransitionTimes', ...
                   sprintf('copied over: start: %s; end: %s', ...
                   nandatestr(datenum(goodStarts(iGoodStart))),...
                   nandatestr(datenum(goodEnds(iGoodEnd)))));
               L.debug('checkTransitionTimes', ...
                   sprintf('copied over: startFlag: %u; endFlag: %u', ...
                   goodStartFlags(iTempStart), goodEndFlags(iTempStart)));

               % update iterators
               iGoodStart = iGoodStart+1;
               iGoodEnd = iGoodEnd+1;

           elseif numEndsFound == 0
              L.debug('checkTransitionTimes', ...
                  'we have no matching end'); 

              % if it's the last one, copy in last ts as end
              if (iTempStart == numTempStarts)
%                   datestr(tempStarts(iTempStart))
                  
                  % <------- COPY IN DATA 2/3 -------------------------->
                  goodStarts(iGoodStart) = tempStarts(iTempStart);
                  goodStartFlags(iGoodStart) = tempStartFlags(iTempStart);

%                   datestr(goodStarts(iGoodStart))
%                
%                   datestr(timestamps(end))
                  goodEnds(iGoodEnd) = timestamps(end);
                  % 11/17 -- here the last timestamp shoudl be one because
                  % we made it up?
%                   datestr(goodEnds(iGoodEnd))
                  goodEndFlags(iGoodEnd) = 1;   % interpolated?

                   L.debug('checkTransitionTimes', ...
                       sprintf('copied over: start: %s; end: %s', ...
                       nandatestr(datenum(goodStarts(iGoodStart))),...
                       nandatestr(datenum(goodEnds(iGoodEnd)))));
                   L.debug('checkTransitionTimes', ...
                       sprintf('copied over: startFlag: %u; endFlag: %u', ...
                       goodStartFlags(iTempStart), goodEndFlags(iTempStart)));
                  
                  % update iterators:
                  iGoodStart = iGoodStart+1;
                  iGoodEnd = iGoodEnd+1;

              else
                  L.debug('checkTransitionTimes', ...
                      'We have a start with no matching end');
                  
                   goodStarts(iGoodStart) = tempStarts(iTempStart);
                   goodStartFlags(iGoodStart) = 3;

                   goodEnds(iGoodEnd) = NaN;
                   goodEndFlags(iGoodEnd) = 3; % interpolated?
                   
                    L.debug('checkTransitionTimes', ...
                       sprintf('copied over: start: %s; end: %s', ...
                       nandatestr(datenum(goodStarts(iGoodStart))),...
                       nandatestr(datenum(goodEnds(iGoodEnd)))));
                    L.debug('checkTransitionTimes', ...
                       sprintf('copied over: startFlag: %u; endFlag: %u', ...
                       goodStartFlags(iTempStart), goodEndFlags(iTempStart)));
                   
                   
                   % update iterators:
                  iGoodStart = iGoodStart+1;
                  iGoodEnd = iGoodEnd+1;                 

              end
           else %numendsfound > 1
                       
               L.debug('checkTransitionTimes', ...
                   'we have more than one match....................');
               % use the first one
               % <------- COPY IN DATA 3/3 -------------------------->
               % first, copy in start time
%                datestr(tempStarts(iTempStart))
               goodStarts(iGoodStart) = tempStarts(iTempStart);
               goodStartFlags(iGoodStart) = tempStartFlags(iTempStart);

               useThisOneIndex = find(possibleEndsIndex == 1, 1, 'first'); 
%                datestr(rawEndTimes(useThisOneIndex))

               goodEnds(iGoodEnd) = rawEndTimes(useThisOneIndex);
               goodEndFlags(iGoodEnd) = tempEndFlags(useThisOneIndex);

               % what do we do with other ones found?  label as 3
               % how many do we have?
               numEndsFound = sum(possibleEndsIndex);
               numExtraEnds = numEndsFound - 1;
               
               % get indices
               extraEndsIndex = find(possibleEndsIndex == 1);
               % get rid of first one
               extraEndsIndex(1) = [];
               % left with what's leftover
                % need to put in nan start, flag as 3, copy over end, flag
                % as 3
               for iExtraEnds = 1:numExtraEnds
                   % 7 + 1
                   goodStarts(iGoodStart + 1) = NaN;
                   
                   % THIS WOULD BE WHERE WE WOULD ADD BITWISE FLAGGING --
                   % WOULD SAY: WHY SUSPECT?
                   goodStartFlags(iGoodStart + 1) = 3;  % suspect
                   goodEnds(iGoodEnd + 1) = rawEndTimes(extraEndsIndex(iExtraEnds));
                   goodEndFlags(iGoodEnd + 1) = 3; % suspect
                   % update interators: plus 2, because we added data to +1
                   iGoodStart = iGoodStart + 1;
                   iGoodEnd = iGoodEnd + 1;
               end   % for loop
               
               iGoodStart = iGoodStart + 1;
               iGoodEnd = iGoodEnd + 1;

           end   % numends found == 1
       else  %is nan
           L.error('checkTransitionTimes', 'if tempStart is nan, why are we using?');
       end
   end   % for loop iStarts
                   
   %  assign out variables
   goodStartsOut = goodStarts;
   goodEndsOut = goodEnds;
   startFlagsOut = goodStartFlags;
   endFlagsOut = goodEndFlags;
             
end   % end function