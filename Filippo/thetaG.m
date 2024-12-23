function [tethaG] = thetaG(year, month, day, hour, minute, second, w)
%THETAG Summary of this function goes here
%   Detailed explanation goes here

w = rad2deg(w) * 3600;
ut = hour + minute/60 + second/3600;
j0 = J0(year, month, day);
j = (j0 - 2451545)/36525;
g0 = 100.4606184 + 36000.77004*j + 0.000387933*j^2- 2.583e-8*j^3;
g0 = wrapTo360(g0);
gst = g0 + w * ut;
tethaG = deg2rad(wrapTo360(gst));
end

