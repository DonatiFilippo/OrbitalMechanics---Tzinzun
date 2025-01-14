function a = a3B_M(t, r, mu_M, date)

% a3B_M - Perturbing acceleration due to the gravitational
% influence of the Moon.
% 
% PROTOTYPE: 
%   a = a3B_M(t, r, mu_M, date)
%
% DESCRIPTION:
%   It gives the perturbing acceleration due to gravitational influence of
%   the Moon acting on a body at distance |r| from Earth, at a given epoch, 
%   expressed in Geocentric Equatorial Reference Frame.
%   This frame {x,y,z} is characterised by:
%       x-axis: on the equatorial plane, along the direction of the gamma
%           point
%       z-axis: direction of the north pole
%       y-axis: on the equatorial plane, completes the reference frame
%
% INPUT: 
%   t [1x1]       Time of evaluation                            [s]
%   r [3x1]       Position of the body at time t in Cartesian 
%                 coordinates                                   [km]
%   mu_M [1x1]    Gravitational parameter of the Moon           [km^3/s^2]   
%	date[1x6]     Date in the Gregorian calendar, as a 6-element vector 
%                 [year, month, day, hour, minute, second]
%
% OUTPUT:
%   a [3x1]       Perturbing acceleration                        [km/s^2]
%
% CONTRIBUTORS:
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo 
%   Domenichelli Eleonora
%
%-------------------------------------------------------------------------

% Conversion of date in Gregorian calendar to date in MJD 2000
MJD2000 = date2mjd2000(date) + t/(60*60*24);

% Computation of Moon's position vector in Cartesian coordinates, expressed
% in Geocentric Equatorial reference frame
[rE_M, ~] = ephMoon(MJD2000);

% Position vector of the body with respect to Moon in Cartesian coordinates
rSC_M = rE_M' - r;

% Perturbing acceleration acting on the bpdy
a = mu_M * ( rSC_M / (norm(rSC_M)^3) - rE_M' / (norm(rE_M')^3) );
end