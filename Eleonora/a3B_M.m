function a = a3B_M(t, r, mu_M, date)

% a3B_M - Perturbing acceleration due to the gravitational
% influence of the Moon.
% 
% PROTOTYPE: 
%   a = a3B_M(t, r, mu_M, date)
%
% DESCRIPTION:
%   It gives the perturbing acceleration due to gravitational influence of
%   the Moon on a at a given epoch in Geocentric Equatorial Reference Frame.
%   This frame {x,y,z} is characterised by:
%       x-axis: on the equatorial plane, along the direction of the gamma
%           point
%       z-axis: direction of the north pole
%       y-axis: on the equatorial plane, completes the reference frame
%INPUT: 
% r [3xn]         Position of the body over time span tv
%                 (rx, ry, rz)                                  [L]
% thetaG0 [1x1]   Greenwich sidereal time at 0 hours UT         []
% tv [1xn]        Time span of orbit propagation                [T]
% wE [1x1]        Earth's rotation velocity                     [rad/T]
% c [1x1]         Parameter to choose GT color                  [-]
%
% OUTPUT:
% alpha [1xn]     Right ascension in ECI in time span tv        [rad]
% delta [1xn]     Declination in ECI in time span tv            [rad]
% lon [1xn]       Longitude with respect to rotating Earth      [deg]
% lat [1xn]       Latitude with respect to rotating Earth       [deg]
%
% CONTRIBUTORS:
% Eleonora Domenichelli
%
% VERSION:
% 2024/11/17: First version
% 2024/12/22: Second version
%-------------------------------------------------------------------------
%   Detailed explanation goes here

% Position and velocity of the Moon in Geocentric equatorial reference
% frame
mjd2000 = date2mjd2000(date);
MJD2000 = mjd2000 + t/(60*60*24);
[rE_M, ~] = ephMoon(MJD2000);

% Position of the spacecraft wrt Moon in Geocentric equatorial reference
% frame
rSC_M = rE_M' - r;

a = mu_M * ( rSC_M / (norm(rSC_M)^3) - rE_M' / (norm(rE_M')^3) );
end