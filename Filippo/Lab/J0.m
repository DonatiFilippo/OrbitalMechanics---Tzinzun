function j0 = J0(year, month, day)
%J0 Summary of this function goes here
%   Detailed explanation goes here
j0 = 367*year - fix(7*(year + fix((month + 9)/12))/4) ...
    + fix(275*month/9) + day + 1721013.5;
end

