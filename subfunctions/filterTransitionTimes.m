function [ goodStartsOut, goodEndsOut ] = ...
    filterTransitionTimes( rawStartsIn, rawEndsIn, durationIn, durToleranceIn, ...
    frequencyIn, freqToleranceIn ) 

% for now lets' jsut start with fresh flags

%FILTERTRANSTIONTIMES - This function checks a set of start and end
%transition times against some rules:
%
% The output flags use 1 == good, 3 == suspect.  

% Any original flags are
% copied over, but if a transition is bad, that flag will be set to
% suspect.
% 
%
% Syntax:  [ goodStartsOut, goodEndsOut, startFlagsOut, endFlagsOut ] = ...
%    filterTransitionTimes( rawStartsIn, rawEndsIn, startFlagsIn, endFlagsIn, timestampsIn )
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
%    filterTransitionTimes( starts, ends, startflags, endflags, timestamps ); 
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
% Apr 2016; Last revision: 03-Apr-16

%------------- BEGIN CODE --------------

L = log4m.getLogger();

% assign input variables
starts = rawStartsIn;
ends = rawEndsIn;
duration = durationIn;   %10 min for FSW; 50 for TSW
frequency = frequencyIn; %60 min for whole cycle
freqTolerance = freqToleranceIn; %
durTolerance = durToleranceIn;

% these are incoming flags that might have been set for a different reason
% startFlags = startFlagsIn;
% endFlags = endFlagsIn;

durationUpperMargin = datenum(0,0,0,0, (duration + durTolerance*duration), 0);
durationLowerMargin = datenum(0,0,0,0, (duration - durTolerance*duration), 0);
frequencyUpperMargin = datenum(0,0,0,0, (frequency) + freqTolerance*frequency, 0);
frequencyLowerMargin = datenum(0,0,0,0, (frequency) - freqTolerance*frequency, 0);

% just for dev:
startFlags = ones(size(starts));
endFlags = ones(size(ends));

jEnds = 1;
lastGoodEnd = 1;
for iStarts = 1:length(starts)
    sprintf('Trying %s', datestr(starts(iStarts)))
    while jEnds <= length(ends)
        periodLength = ends(jEnds) - starts(iStarts);
        % if this period is outside margins for period duration
        if (periodLength < durationLowerMargin ) %|| periodLength > durationUpperMargin)
%         if (periodLength < durationLowerMargin) % if it's too short
            sprintf('too short: %s'), datestr(ends(jEnds))
            % not a good end, mark as 3.
           
            endFlags(jEnds) = 3;
            startFlags(iStarts) = 3;
%         jEnds = jEnds+1;
            break;
        elseif (periodLength > durationUpperMargin)&& jEnds == length(ends)
            sprintf('too long and LAST ONE: %s'), datestr(ends(jEnds))
            % not a good end, mark as 3.
           % let's not flag it for now
            endFlags(jEnds) = 3;
            startFlags(iStarts) = 3;
            break;
        else
            % should be good?
            sprintf('good: %s'), datestr(ends(jEnds))
            endFlags(jEnds) = 1; %won't revisit
            startFlags(iStarts) = 1;
            lastGoodEnd = jEnds;
            break;
        end;
        % if we get here, we found no good end for this one, mark start as
        % bad
%         startFlags(iStarts) = 3;
%         jEnds = jEnds+1;
    end;
    jEnds = jEnds +1;
    disp('-----------');
end;
    
goodStartsOut = starts(startFlags==1);
goodEndsOut = ends(endFlags==1);


    