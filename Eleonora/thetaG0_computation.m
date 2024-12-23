function thetaG0 = thetaG0_computation(J0, UT, wE)
%UNTITLED4 Summary of this function goes here
% wE is in rad/s but needs to be in deg/h

wE = rad2deg(wE) * 3600;
T0 = (J0 - 2451545)/36525;
thetaG_t0 = wrapTo360(100.4606184 + 36000.77004*T0 + 0.000387933*(T0^2) - 2.583*(10^-8)*(T0^3));
thetaG0 = wrapTo360(thetaG_t0 + wE*UT);
thetaG0 = deg2rad(thetaG0);
end