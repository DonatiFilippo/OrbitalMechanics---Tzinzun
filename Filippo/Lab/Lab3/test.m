clc
close all
clear all

mu_E = astroConstants(13);
om_E = deg2rad(15.4) * (1/3600);
[r0, v0] = parOrb2rv(7091.185, 0, deg2rad(30), 0, deg2rad(40), deg2rad(0));
y0 = [r0; v0];

a = 1/(2/norm(r0) - dot(v0, v0)/mu_E);
T = 2*pi*sqrt(a^3/mu_E);
tspan = linspace(0, 10*T, 300000);

% %% Plot
% close all
% 
% img = imread('EarthTexture.jpg');
% image(img);
% hold on;
% plot(v*100);

%% test

[alpha, phi, lon, lat] = groundTrack2(y0, tspan, 0, mu_E, om_E);
