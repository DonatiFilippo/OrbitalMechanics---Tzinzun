function a = aJ2(r, mu, R, J2)

% aJ2 - Perturbing acceleration due to second zonal harmonic J2.
%
% PROTOTYPE:
%   a = aJ2(r, mu, R, J2)
%
% DESCRIPTION:
%   The function evaluates the perturbing acceleration due to second zonal 
%   harmonic J2 on a body at distance |r| from teh primary, expressed in
%   Geocentric Equatorial reference frame.
%   This frame {x,y,z} is characterised by:
%       x-axis: on the equatorial plane, along the direction of the gamma
%           point
%       z-axis: direction of the north pole
%       y-axis: on the equatorial plane, completes the reference frame
%
% INPUT:
%   r[3x1]   Position of the body in Cartesian coordinates       [km]
%   mu[1x1]  Gravitational parameter of the primary              [km^3/s^2]
%   R[1x1]   Radius of the primary                               [L]
%   J2[1x1]  Gravitational harmonic coefficient of the primary   [-]
%
% OUTPUT: 
%    a[3X1]  Perturbing acceleration due to J2                   [km/s^2]
%
% CONTRIBUTORS:
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo 
%   Domenichelli Eleonora
%
%-------------------------------------------------------------------------

% Distance from the primary
r_norm = norm(r);

% Perturbing acceleration
a = 3/2*( J2*mu*R^2 )/r_norm^4 * [r(1)/r_norm * (5*(r(3)/r_norm)^2-1);...
                                 r(2)/r_norm * (5*(r(3)/r_norm)^2-1);...
                                 r(3)/r_norm * (5*(r(3)/r_norm)^2-3)];

end