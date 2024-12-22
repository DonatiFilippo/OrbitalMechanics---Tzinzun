function [r_vect,v_vect] = kep2cart (a,e,i,OM,om,th,mu)
%
% Converter: Keplerian Parameters ---> Cartesian Coordinates & Velocities
%
%DESCRIPTION:
%This code provides a conversion from the Keplerian Parameters of the orbit
%to cartesian coordinates vector (r_vect) and velocities (v_vect).
%
%--------------------------------------------------------------------------
% INPUTS:
%   a          [1x1]       Semi-Major Axis           [km]
%   e          [1x1]       Eccentricity              [-]
%   i          [1x1]       Inclination               [rad]
%   OM         [1x1]       RAAN                      [rad]
%   om         [1x1]       Argument of Periapsis     [rad]
%   th         [1x1]       True Anomaly              [rad]
%   mu         [1x1]       Planetary Constant        [km^3][sec-2]
%--------------------------------------------------------------------------
% OUTPUTS:
%   r_vect     [3x1]       Position Vector           [km]
%   v_vect     [3x1]       Velocity Vector           [km/s]
%--------------------------------------------------------------------------
%
%% Conversion Routine
%Semi-Latus Rectum
p = a*(1-e.^2);
%Radial Distance
r=p/(1+e*cos(th));
%Position and Velocity Vectors in Perifocal Coordinates
r_pf = r*[cos(th);sin(th);0];
v_pf = sqrt(mu/p)*[-sin(th);e+cos(th);0];
%Rotation Matrix
R = [cos(om)*cos(OM)-sin(om)*cos(i)*sin(OM) -sin(om)*cos(OM)-cos(om)*cos(i)*sin(OM) sin(i)*sin(OM);
     cos(om)*sin(OM)+sin(om)*cos(i)*cos(OM) -sin(om)*sin(OM)+cos(om)*cos(i)*cos(OM) -sin(i)*cos(OM);
     sin(om)*sin(i)                         cos(om)*sin(i)                          cos(i)];
%Position and Velocity Vectors in Cartesian Coordinates
r_vect = R*r_pf;
v_vect = R*v_pf;
end