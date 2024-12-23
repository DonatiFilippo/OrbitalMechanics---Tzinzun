function [J0, UT] = J0_computation(jd, date)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
UT = date(4) + date(5)/60 + date(6)/3600;
J0 = jd - UT/24;
end