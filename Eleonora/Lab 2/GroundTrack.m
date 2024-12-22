function [alpha, delta, lon, lat] = GroundTrack(y0, thetaG0, tv, wE, mu)
% Computation of latitude and longitude over time interval tv
%
% PROTOTYPE: [alpha, delta, lon, lat] = GroundTrack(y0, thetaG0, t, wE, mu)
%
% INPUT: 
% y0 [6x1]
% OUTPUT:
%
% CONTRIBUTORS:
%  Eleonora Domenichelli
%
% VERSION:
%  2024/11/17: First version
%
%-------------------------------------------------------------------------

% First we need to propagate our orbit using function 2bp implemented in
% last lab

% Options for the ODE solver
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

% Perform the integration
[~, Yout] = ode113 (@(t,y) ode_2bp(t,y,mu), tv, y0, options);
Yout = Yout'; % [6xl] row vector

% Compute right ascension and declination in time
l = length(tv);
delta = zeros(1, l);
alpha = zeros (1, l);
lon = zeros (1, l);

for i = 1:l
    r = norm(Yout(1:3, i));
    delta(i) = asin(Yout(3,i)/r);
    alpha(i) = atan2(Yout(2,i), Yout(1, i));

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

