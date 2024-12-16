function [alpha,delta, lon, lat] = groundTrack2(y0, tv, theta, mu, om)
%GROUNDTRACK Summary of this function goes here
%   Detailed explanation goes here


%% Propagator
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[tu, yu] = ode113(@(t, y) ode_2bp(t, y, mu), tv, y0, options);
rr = yu(:, 1:3);
r = vecnorm(rr, 2, 2);

%% Conversion to RA and Declination
delta = asin(rr(:, 3)./r);
alpha = atan2(rr(:,2), rr(:,1));

%% Conversion to Latitude and Longitude
g = theta + om .* tu;
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


