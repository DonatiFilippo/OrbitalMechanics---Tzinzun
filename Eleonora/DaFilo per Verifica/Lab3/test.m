clc
close all
clear all

mu_E = astroConstants(13);
om_E = deg2rad(15.04) * (1/3600);
a = 42166.167;
[r0, v0] = kep2car(a, 0, deg2rad(0), deg2rad(0), deg2rad(0), deg2rad(20),mu_E);
y0 = [r0; v0];

T = 2*pi*sqrt(a^3/mu_E);
tspan = linspace(0, 30*T, 3000000);

% %% Plot
% close all
% 
% img = imread('EarthTexture.jpg');
% image(img);
% hold on;
% plot(v*100);

%% test

[alpha, phi, lon, lat] = groundTrack2(y0, tspan, 0, mu_E, om_E);
