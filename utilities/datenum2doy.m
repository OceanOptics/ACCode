function doy = datenum2doy(dn)
dv = datevec(dn);
dv(:, 2:end) = 0;
doy = dn - datenum(dv);
end