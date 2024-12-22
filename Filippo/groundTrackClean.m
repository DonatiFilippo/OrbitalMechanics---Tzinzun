function [alpha,delta, lon, lat] = groundTrackClean(s, t, theta, om0)
%groundTrack Evaluate and plots the ground track of an object
%   
% PROTOTYPE: [alpha,delta, lon, lat] = groundTrackClean(s, t, theta, om0)
%
% INPUT:
%   s[nx6]      Cartesian state of the object (rx, ry, rz, vx, vy, vz) [m, m/s]
%   t[nx1]      Time Vector [s]
%
% OUTPUT:
%
%
% CONTRIBUTORS:
%   Azevedo Da Silva Esteban
%   Domenichelli Eleonora
%   Donati Filippo
%   Gavidia Pantoja Maria Paulina
%
% VERSIONS
%   18-12-2024: Initial Commit

%% Variable extraction

rr = s(:, 1:3);
r = vecnorm(rr, 2, 2);

%% Conversion to RA and Declination
delta = asin(rr(:, 3)./r);
alpha = atan2(rr(:,2), rr(:,1));

%% Conversion to Latitude and Longitude
g = theta + om0 .* t;
lat = delta;
lon = alpha - g;
lon = wrapToPi(lon);

%% Plot
figure
lon_deg = rad2deg(lon);
lat_deg = rad2deg(lat);
hold on;
img = imread('EarthTexture.jpg');
image([-180,+180], [+90, -90],img);
plot(lon_deg, lat_deg, 'LineStyle','none','Marker','.', 'Color','0 1 0','DisplayName', 'Ground Track');
plot(lon_deg(1,1), lat_deg(1,1), 'LineStyle','none','Marker','o','DisplayName', 'Start'); % Start Point
plot(lon_deg(end,1), lat_deg(end,1), 'LineStyle','none','Marker','square', 'DisplayName', 'End'); % End Point
xlabel('Longitude [deg]');
xlim([-180 +180]);
ylim([-90 +90]);
ylabel('Latitude [deg]');
legend;
grid on;
end

