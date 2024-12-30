function [rr, vv] = parorb2rv (a, e, i, OM, om, theta, mu)

% parorb2rv - Converts Keplerian elements to Cartesian state vectors.
%
% PROTOTYPE:
%   [rr, vv] = parorb2rv(a, e, i, OM, om, theta, mu)
%
% DESCRIPTION:
%   The function converts Keplerian elements to Cartesian state vector.
%   Inputs can be only scalar.
%
% INPUT:
%   a [1x1]       Semi-major axis                             [km]
%   e [1x1]       Eccentricity                                [-]
%   i [1x1]       Inclination                                 [rad]
%   OM [1x1]      Right ascension of ascending node           [rad]
%   om [1x1]      Pericenter anomaly                          [rad]
%   theta [1x1]   True anomaly                                [rad]
%   mu [1x1]      Gravitational parameter of the primary      [km^3/s^2]
%
% OUTPUT:
%   rr [3x1]      Position vector in Cartesian coordinates    [km]
%   vv [3x1]      Velocity vector in Cartesian coordinates    [km/s]
%
% CONTRIBUTORS:
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo 
%   Domenichelli Eleonora
%
%-------------------------------------------------------------------------

% Semi-latus rectum
p=a*(1-(e^2));

% Position and velocity vectors expressed in perifocal frame
rr_PF=(p/(1+e*cos(theta)))*[cos(theta); sin(theta); 0];
vv_PF=sqrt(mu/p)*[-sin(theta); e+cos(theta); 0];

% Rotation matrix from perifocal to ECI definition
R_OM=[cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R_i=[1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R_om=[cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];

T_PF_GE=(R_om*R_i*R_OM)';

% Rotation of position and velocity vectors
rr=T_PF_GE*rr_PF;
vv=T_PF_GE*vv_PF;
end