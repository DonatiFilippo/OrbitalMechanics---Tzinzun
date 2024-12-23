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
% Compute right ascension and declination in time
delta = zeros(1, n);
alpha = zeros (1, n);
lon = zeros (1, n);

for i = 1:n
    r_norm = norm(r(:, i));
    delta(i) = asin(r(3,i)/r_norm);
    alpha(i) = atan2(r(2,i), r(1, i));

    % Longitude
    thetaG = thetaG0 + wE*tv(i);
    lon(i) = (alpha(i) - thetaG);

    if ((lon(i) < -pi) || (lon(i) > pi))
        lon(i) = mod(lon(i)+pi, 2*pi) - pi; % to put it in the interval -pi, pi
    end

end

% Compute latitude 
lat = delta;

% Rad to deg conversion for plot and results readability
lon = rad2deg(lon);
lat = rad2deg(lat);

% Ground track plot
S = imread("EarthTexture.jpg");

figure
image([-180, 180], [90, -90], S);
hold on;
grid minor;

% Green ground track is the non repeated one, red is the repeated one
if c == 1 
    plot(lon, lat, 'g', 'LineStyle','none','Marker','.');
else 
    plot(lon, lat, 'r', 'LineStyle','none','Marker','.');
end

plot(lon(1,1), lat(1,1), 'ko', 'LineWidth', 2)
plot(lon(1,end), lat(1,end), 'ks', 'LineWidth', 2)
set(gca, 'YDir', 'normal')

xlabel("Longitude [deg]");
ylabel("Latitude [deg]");
legend("Ground track", "Start", "End");
end

% TOGLIERE C PERCHE' è INUTILE