function a3B_M = a3B_M(t, r, mu_M, date)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Position and velocity of the Moon in Geocentric equatorial reference
% frame
mjd2000 = date2mjd2000(date);
MJD2000 = mjd2000 + t/(60*60*24);
[rE_M, ~] = ephMoon(MJD2000);

% Position of the spacecraft wrt Moon in Geocentric equatorial reference
% frame
rSC_M = rE_M' - r;

a3B_M = mu_M * ( rSC_M / (norm(rSC_M)^3) - rE_M' / (norm(rE_M')^3) );
end