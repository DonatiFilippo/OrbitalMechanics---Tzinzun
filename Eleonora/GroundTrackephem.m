function [alpha, delta, lon, lat] = GroundTrackephem(y0, thetaG0, tv, wE)
% Computation of latitude and longitude using ephemerides
%
% PROTOTYPE: [alpha, delta, lon, lat] = GroundTrackephem(y0, thetaG0, t, wE, mu)
%
% INPUT: 
% y0 [6xl]
% OUTPUT:
%
% CONTRIBUTORS:
%  Eleonora Domenichelli
%
% VERSION:
%  2024/11/17: First version
%
%-------------------------------------------------------------------------

% With ephemerides I don't need to propagate my orbit, already done

% Compute right ascension and declination in time
l = length(y0);
delta = zeros(1, l);
alpha = zeros (1, l);
lon = zeros (1, l);

for i = 1:l
    r = norm(y0(1:3, i));
    delta(i) = asin(y0(3,i)/r);
    alpha(i) = atan2(y0(2,i), y0(1, i));

    % Longitude
    thetaG = thetaG0 + wE*tv(i);
    lon(i) = (alpha(i) - thetaG);

    if ((lon(i) < -pi) || (lon(i) > pi))
        lon(i) = mod(lon(i)+pi, 2*pi) - pi; % to put it in the interval -pi, pi
    end

end

% Compute latitude 
lat = delta;
end
