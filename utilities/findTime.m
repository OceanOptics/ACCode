% find a row based on a timestamp

function [ind] = findTime(bins, year, mon, day, hour, min, sec)

    find_time = datenum(year, mon, day, hour, min, sec)
    [ind] = datefind(find_time, bins, 0)
end
