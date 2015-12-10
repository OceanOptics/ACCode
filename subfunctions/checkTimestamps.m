function boolean_result = checkTimestamps(ts1, ts2)

    boolean_result = true;
    if size(ts1) ~= size(ts2)
%         disp('size not same')
        boolean_result = false;
    else
%         disp('size same')
        % size is the same
        if sum(datevec(ts1(1)) == datevec(ts2(1))) == 6
%             disp('same first ts')
        else
%             disp('diff first ts')
            boolean_result = false;
        end;
        if sum(datevec(ts1(end)) == datevec(ts2(end))) == 6
%             disp('same end ts')
        else
%             disp('diff last ts')
            boolean_result = false;
        end;
    end;