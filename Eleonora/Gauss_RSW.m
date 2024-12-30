function dkep = Gauss_RSW(t, kep, mu, R, J2, muM, date)

% Gauss_RSW - Gauss planetary equation with perturbing acceleration in RSW 
% (radial-transversal-out of plane) reference frame
%
% PROTOTYPE:
%   dkep = Gauss_RSW(t, kep, mu, R, J2, muM, date)
% 
% DESCRIPTION:
%   The function gives Gauss planetary equations (set of differential
%   equations that describe the time evolution of the orbital elements of a
%   celestial body under the influence of perturbing accelerations) in the
%   formulation with perturbing acceleration in RSW (radial-transversal-out
%   of plane) reference frame. The perturbing accelerations are computed
%   directly inside this function using the auxiliary functions a3B_M and
%   aJ2.
%   {r, w, s} is characterized by:
%       r-axis: on the orbit plane, along the direction of the radius
%       w-axis: out of plane, along the direction of the angular momentum
%       of the orbit
%       s-axis: on the orbit plane, completes the reference frame
%
% INPUT:
%   t [1x1]       Time instant for acceleration evaluation       [s]
%   kep [1x6]     Keplerian elements vector
%                 (a, e, i, OM, om, theta)                       [km, -, rad] 
%   mu [1x1]      Gravitational parameter of the primary         [km^3/s^2]
%   R [1x1]       Radius of the primary                          [km]
%   J2 [1x1]      Gravitational harmonic coefficient of
%                 the primary                                    [-]
%   muM [1x1]     Gravitational parameter of the Moon            [km^3/s^2]
%
%   date [1x6]    Date in the Gregorian calendar, as a 6-elements vector
%                 [year, month, day, hour, minute, second]. For dates 
%                 before 1582, the resulting date components are valid 
%                 only in the Gregorian proleptic calendar. This is based
%                 on the Gregorian calendar but extended to cover dates 
%                 before its introduction. Date must be after 12:00 noon,
%                 24 November -4713.
%
% OUTPUT:
%   dkep [6x1]    Gauss' planetary equations                      [km/s, s^-1, rad/s]
%
% CONTRIBUTORS:
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo 
%   Domenichelli Eleonora
%
%-------------------------------------------------------------------------

% Orbital parameters extraction
a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
theta = kep(6);

% Conversion of Keplerian parameters to radius and velocity vectors
% in Cartesian coordinates
[rp, vp] =  parorb2rv (a, e, i, OM, om, theta, mu);

% Evaluate perturbing acceleration, expressed in Geocentric equatorial
% reference frame
ap_ECI = aJ2(rp, mu, R, J2) + a3B_M(t, rp, muM, date);

% Acceleration's rotation in RSW reference frame
ap = TNH2RSW(kep, mu) * GE2TNH(rp, vp) * ap_ECI;

% Computation of useful values for equations readability
p = a*(1-e^2);
r = p / ( 1+e*cos(theta));
n = sqrt(mu/(a^3));
b = a*sqrt(1-e^2);
h = n * a * b;

% Keplerian parameters evolution
dkep = [2*(a^2)/h * (e*sin(theta)*ap(1,1) + p/r * ap(2,1));
        1/h * (p*sin(theta)*ap(1,1) + ((p+r)*cos(theta) + r*e)*ap(2,1));
        r*cos(theta+om)/h * ap(3,1);
        r*sin(theta+om)/(h*sin(i)) * ap(3,1);
        1/(h*e) * (-p*cos(theta)*ap(1,1) + (p+r)*sin(theta)*ap(2,1)) - (r*sin(theta+om)*cos(i))/(h*sin(i)) * ap(3,1);
        h/(r^2) + 1/(e*h) * (p*cos(theta)*ap(1,1) - (p+r)*sin(theta)*ap(2,1))];
end