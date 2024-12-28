function A_TNH = GE2TNH(y)

% GE2TNH - Rotation matrix from Geocentric Equatorial reference frame to
% TNH (tangential-normal-out of plane) reference frame
% 
% PROTOTYPE: 
%   A_TNH = car2TNH(y)
%
% DESCRIPTION:
%   The function gives the rotation matrix from Geocentric Equatorial
%   reference frame to TNH reference frame, taking as input the state
%   (position and velocity) of the body, at the time istant we want to
%   operate the reference frame change.
%   The first frame {x,y,z} is characterised by:
%       x-axis: on the equatorial plane, along the direction of the gamma
%           point
%       z-axis: direction of the north pole
%       y-axis: on the equatorial plane, completes the reference frame
%   The second frame {t, n, h} is characterised by:
%       t-axis: on the orbit plane, along the direction of the velocity
%       h-axis: out of plane, along the direction of the angular momentum
%       of the orbit
%       n-axis: on the orbit plane, completes the reference frame
%
% INPUT: 
%   y [6x1]       State of the body in Cartesian coordinates
%                 (rx, ry, rz, vx, vy, vz)                       [km, km/s]
%
% OUTPUT:
%   A_TNH [3x3]    Rotation matrix                               [km/s^2]
%
% CONTRIBUTORS:
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo 
%   Domenichelli Eleonora
%
%-------------------------------------------------------------------------

% Variables separation
r = y(1:3,1);
v = y(4:end, 1);

% Tangential unit vector
t_uv = v / norm(v);

% Out of plane unit vector
h = cross(r, v);
h_uv = h / norm(h);

% Normal unit vector
n_uv = cross(h_uv, t_uv);

% Rotation matrix assembly
A_TNH = [t_uv, n_uv, h_uv]';
end