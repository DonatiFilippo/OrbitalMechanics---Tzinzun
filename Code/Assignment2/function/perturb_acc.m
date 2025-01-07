function ap = perturb_acc(t, r, mu, R, J2, mu_3B)
% Evaluate the perturbing acceleration due to second zonal harmonic J2,
% expressed in Geocentric equatorial reference frame
%
% PROTOTYPE:
%  a = perturb_acc(y, mu, R, J2)
%
% INPUT:
% r[3x1]     Position of the body in ECI (rx, ry, rz)            [L, L/T]
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

if nargin == 5 %only computation of J2 to TEST
% Distance from the primary
r_norm = norm(r);

% Perturbing acceleration due to J2
aJ2 = 3/2*( J2*mu*R^2 )/r_norm^4 * [r(1)/r_norm * (5*(r(3)/r_norm)^2-1);...
                                 r(2)/r_norm * (5*(r(3)/r_norm)^2-1);...
                                 r(3)/r_norm * (5*(r(3)/r_norm)^2-3)];

ap = aJ2;

else
% Distance from the primary
r_norm = norm(r);

% Perturbing acceleration due to J2
aJ2 = 3/2*( J2*mu*R^2 )/r_norm^4 * [r(1)/r_norm * (5*(r(3)/r_norm)^2-1);...
                                 r(2)/r_norm * (5*(r(3)/r_norm)^2-1);...
                                 r(3)/r_norm * (5*(r(3)/r_norm)^2-3)];
% Position and velocity of the Moon in Geocentric equatorial reference
% frame
[rE_M, ~] = ephMoon(MJD2000);

% Position of the spacecraft wrt Moon in Geocentric equatorial reference
% frame
rSC_M = rE_M - r;

% Perturbing acceleration due to moon
a3B_M = mu_3B * ( rSC_M / (norm(rSC_M)^3) - rE_M / (norm(rE_M)^3) );

ap = aJ2 + a3B_M;
end

end
