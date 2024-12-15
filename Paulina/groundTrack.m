function [alpha, delta, lon, lat] = groundTrack(y0, thetaG0, t, mu, omega_E, t0)

% groundTrack.m - Computes the ground track of a satellite in orbit.
% 
% Inputs:
%   r_0       - Initial position vector in Cartesian coordinates [x, y, z] (m)
%   theta_G0  - Longitude of Greenwich meridian at initial time (rad)
%   times     - Vector of times at which the ground track will be computed (s)
%   mu        - Gravitational parameter of Earth (m^3/s^2)
%   omega_E   - Earth's angular velocity (rad/s)
%   t0        - Initial time (s)
%
% Outputs:
%   alpha     - Right ascension in Earth-centered inertial frame (rad)
%   delta     - Declination in Earth-centered inertial frame (rad)
%   lon       - Longitude with respect to rotating Earth (rad)
%   lat       - Latitude with respect to rotating Earth (rad)
%
% Author
%  Maria Paulina Pantoja Gavidia
% ------------------------------------------------------------------------

% Extract initial position and velocity
r0 = y0(1:3);
v0 = y0(4:6);

% Propagate orbit using ODE45 (two-body dynamics)
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-9);
[~, y] = ode45(@(t, y) ode_2bp(t, y, mu), t, y0, options);

% Extract position vectors at each time
r_vectors = y(:, 1:3);

% Initialize outputs
alpha = zeros(length(t), 1);
delta = zeros(length(t), 1);
lon = zeros(length(t), 1);
lat = zeros(length(t), 1);

% Compute alpha, delta, longitude, and latitude
for i = 1:length(t)
    r = r_vectors(i, :);
    r_norm = norm(r);
    
    % Declination
    delta(i) = asin(r(3) / r_norm);
    
    % Right ascension
    alpha(i) = atan2(r(2), r(1));
    
    % Greenwich Sidereal Time
    thetaG = thetaG0 + omega_E * (t(i) - t0);
    
    % Longitude (adjusted for Earth's rotation)
    lon(i) = wrapToPi(alpha(i) - thetaG);
    
    % Latitude
    lat(i) = delta(i);
end

end