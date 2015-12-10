function [start_time_index end_time_index] = findDateIndices(cu_time, start_time, end_time)
    
    %findDateIndices
    %Finds the indices of a time array corresponding to the closest times to
    %start_time and end_time.
    
    start_time_index = min(find( min(abs(cu_time - start_time)) == abs(cu_time - start_time) , 1 ));
    end_time_index = max( find( min(abs(cu_time - end_time)) == abs(cu_time - end_time) , 1, 'last' ));

end

% (cu_time - start_time) will give back an index of numbers of distance
% between each timestam and the specified time

% (abs) gives the absolute value

% min finds the smallest
% 
% so find the number with the minimum abs value of difference
% 
% ,1, 'last' finds the LAST one that matches?

% final min() is to grab the earliest if there is more than one timestamp
% returned.