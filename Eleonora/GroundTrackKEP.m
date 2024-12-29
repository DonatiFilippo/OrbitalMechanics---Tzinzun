function [alpha, delta, lon, lat] = GroundTrackKEP(kep, thetaG0, tv, wE, muE, c)
% Computation of latitude and longitude over time interval tv and ground
% track plot
%
% PROTOTYPE: [alpha, delta, lon, lat] = GroundTrack(y, thetaG0, t, wE, c)
%
% INPUT: 
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

n = length(tv);
r = zeros(3,n);

% Obtain r
for k = 1:n
    r(:,k) = parorb2rv (kep(1,k), kep(2,k), kep(3,k), kep(4,k), kep(5,k), kep(6,k), muE);
end



