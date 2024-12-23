function aJ2 = aJ2(r, mu, R, J2)
% Evaluate the perturbing acceleration due to second zonal harmonic J2,
% expressed in Geocentric equatorial reference frame
%
% PROTOTYPE:
%  a = perturb_acc(y, mu, R, J2)
%
% INPUT:
% y[6x1]     Position of the body in ECI (rx, ry, rz)            [L, L/T]
% mu[1x1]    Gravitational parameter of the primary              [L^3/T^2]
% R[1x1]     Primary's radius                                    [L]
% J2[1x1]    Gravitational harmonic coefficient of the primary   [-]
%
% OUTPUT: 
% a[3X1]     Perturbing acceleration due to J2                    [L/T^2]
%
% CONTRIBUTORS:
%  Eleonora Domenichelli
%
% VERSION:
%  2024/12/22: First version
%
%-------------------------------------------------------------------------

% Distance from the primary
r_norm = norm(r);

% Perturbing acceleration due to J2
aJ2 = 3/2*( J2*mu*R^2 )/r_norm^4 * [r(1)/r_norm * (5*(r(3)/r_norm)^2-1);...
                                 r(2)/r_norm * (5*(r(3)/r_norm)^2-1);...
                                 r(3)/r_norm * (5*(r(3)/r_norm)^2-3)];

end