function [alpha,delta, lon, lat] = groundTrackJ2G(y0, tv, theta, mu, om, J, R, muM, start)
%GROUNDTRACK Summary of this function goes here
%   Detailed explanation goes here


%% Propagator
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[tu, yu] = ode113(@(t, y) gaussFun(t, y, accPert(t, y,mu, J, R, muM,start),mu), tv, y0, options);
l = length(tv);
rr = zeros(l, 3);
for i = 1:l
    rrv = kep2car(yu(i, 1), yu(i, 2), yu(i, 3), yu(i, 4), yu(i, 5), yu(i, 6), mu);
    rr(i,:) = rrv';
end
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

