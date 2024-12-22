%% Lab 2
% Add comments and plot inside groundtrack function, add variable to choose colors
% comment rv2parorb and parorb2rv 
%% Ex 1
clear;
clc;
close all;


r0 = [3108.128, -1040.299, -6090.022];
v0 = [5.743, 8.055, 1.555];
a = 26600;
thetaG0 = 0;
mu = astroConstants(13);
T = 2*pi*sqrt(a^3/mu);
t0 = 0;
tv = linspace(t0, 30*T, 100000);
wE = 15.04 * pi/180*1/3600;
y0 = [r0'; v0'];

[alpha2,delta2, lon2, lat2] = GroundTrack2(y0, tv, thetaG0, mu, wE);
[alpha, delta, lon, lat] = GroundTrack(y0, thetaG0, tv, wE, mu);

% Ground track plot
S = imread("EarthTexture.jpg");
lon = rad2deg(lon);
lat = rad2deg(lat);

figure
image([-180, 180], [90, -90], S);
hold on;
grid minor;
plot(lon, lat, 'g', 'LineStyle','none','Marker','.');
plot(lon(1,1), lat(1,1), 'ro', 'LineWidth', 2)
plot(lon(1,end), lat(1,end), 'rs', 'LineWidth', 2)
set(gca, 'YDir', 'normal')

xlabel("Longitude [deg]");
ylabel("Latitude [deg]");
legend("Ground track", "Start", "End");

% Molniya
% r0 = [3108.128, -1040.299, -6090.022];
% v0 = [5.743, 8.055, 1.555];
% a = 26600;

% General orbit
% r0 = [-4578.219, -801.084, -7929.708];
% v0 = [0.800, -6.037, 1.385];
% a = 8350;

% LEO i=0°
% r0 = [5493.312, 4609.436, 0.000];
% v0 = [-4.792, 5.711, 0.000];
% a = 7171.010;

% LEO i=30°
% r0 = [5493.312, 3991.889, 2304.718];
% v0 = [-4.792,  4.946, 2.856];
% a = 7171.010;

% LEO i=98°
% r0 = [5493.312, -641.510, 4564.578];
% v0 = [-4.792, -0.795, 5.656];
% a = 7171.010;
%% Ex 2
clear;
clc;
close all;

r0 = [3108.128, -1040.299, -6090.022];
v0 = [5.743, 8.055, 1.555];
y0 = [r0'; v0'];
a = 26600;

e = 0.74;
i = deg2rad(63.4);
OM = deg2rad(50);
om = deg2rad(280);
theta = deg2rad(0);
thetaG0 = 0;

mu = astroConstants(13);
T = 2*pi*sqrt(a^3/mu);
t0 = 0;
tv = linspace(t0, 30*T, 10000);
wE = 15.04 * pi/180*1/3600;
k=2;
m=1;

ar = repeatingGroundTrack(m, k, wE, mu);
Tr =  2*pi*sqrt(ar^3/mu);
tvr = linspace(t0, 30*Tr, 10000);
[r0r, v0r] = parorb2rv (ar, e, i, OM, om, theta, mu);
y0r = [r0r'; v0r'];
[alpha, delta, lon, lat] = GroundTrack(y0, thetaG0, tv, wE, mu);

[alphar, deltar, lonr, latr] = GroundTrack(y0r, thetaG0, tvr, wE, mu);

% Ground track plot
S = imread("EarthTexture.jpg");
lon = rad2deg(lon);
lat = rad2deg(lat);
lonr = rad2deg(lonr);
latr = rad2deg(latr);

figure
image([-180, 180], [90, -90], S);
hold on;
grid minor;
plot(lon, lat, 'g', 'LineStyle','none','Marker','.');
plot(lonr, latr, 'r', 'LineStyle','none','Marker','.');
plot(lon(1,1), lat(1,1), 'go', 'LineWidth', 2)
plot(lon(1,end), lat(1,end), 'gs', 'LineWidth', 2)
plot(lonr(1,1), latr(1,1), 'ro', 'LineWidth', 2)
plot(lonr(1,end), latr(1,end), 'rs', 'LineWidth', 2)
set(gca, 'YDir', 'normal')

xlabel("Longitude [deg]");
ylabel("Latitude [deg]");
legend("Ground track", "Repeated GT");

%% Check if parorb2rv and rv2parorb works
clear;
clc;
close all;

mu = astroConstants(13);
a0 = 26600;
e0 = 0.74;
i0 = deg2rad(63.4);
OM0 = deg2rad(50);
om0 = deg2rad(280);
theta0 = deg2rad(0);

[rr, vv] = parorb2rv (a0, e0, i0, OM0, om0, theta0, mu) %funziona
[a, e, i, OM, om, theta] = rv2parorb(rr, vv, mu)
 % funziona ma nel caso con i=0 e e=0 om e theta sono invertiti

%% Plot ground track from ephemerides
% DON'T WORKS
wE = 15.04 * pi/180*1/3600;
mu = astroConstants(13);
thetaG0 = 0;

tv = [0:300:86400];
y0 = horizonsresults(:,1:3)';
[~, ~, lon, lat] = GroundTrackephem(y0, thetaG0, tv, wE);

S = imread("EarthTexture.jpg");
lon = rad2deg(lon);
lat = rad2deg(lat);

figure
image([-180, 180], [90, -90], S);
hold on;
grid minor;
plot(lon, lat, 'g', 'LineStyle','none','Marker','.');
plot(lon(1,1), lat(1,1), 'ro', 'LineWidth', 2)
plot(lon(1,end), lat(1,end), 'rs', 'LineWidth', 2)
set(gca, 'YDir', 'normal')

xlabel("Longitude [deg]");
ylabel("Latitude [deg]");
legend("Ground track", "Start", "End");

% better in the other way, propagating the time interval, I think thetaG0
% needs a correction
