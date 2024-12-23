%% Main script for Planetary Explorer Mission assignment %%
%
% Group number : 2427
%
%--------------------------------------------------------------------------
% Created and maintained by : 
%
% Azevedo Da Silva Esteban
% Gavidia Pantoja Maria Paulina
% Donati Filippo 
% Domenichelli Eleonora
% 
%--------------------------------------------------------------------------
% LAST UPDATE: 22-12-2024
%

% 2427 / Earth / J2 & Moon / STARTING DATE OF PROPAGATION? 

%% SCRIPT INITIALISATION
clc
clear variables
close all

%% PRESENTATION
fprintf('----------------------------------------------\n');
fprintf('   PLANETARY EXPLORER MISSION ASSIGNMENT      \n');
fprintf('----------------------------------------------\n');
fprintf('\n');

%% Important Constants
mu_E = astroConstants(13);
mu_Moon = astroConstants(20); 
R_E = astroConstants(23);
J2 = astroConstants(9);
w_E = deg2rad(15.04) * (1/3600);

%% ORBIT CHARACTERISATION
a = 23300; % Semi-major axis [km]
e = 0.523; % Eccentricity [-]
i = deg2rad(58.65); % Inclination [rad]
OM = 0; % RAAN [rad]
om = 0; % Argument of pericenter [rad]
theta0 = 0; % True Anomaly [rad]

date = [204, 12, 22, 23, 20, 30];

Earth_3D;
hold on
plotOrbit(a, e, i, OM, om);

k = 13;
m = 3;

%% Ground Track

N = 30000;

T = 2*pi*sqrt(a^3/mu_E);
t1 = linspace(0, T, N);
t2 = linspace(0, 24*60*60, N);
t3 = linspace(0, 30*T, N);

[rr, vv] = kep2car(a, e, i, OM, om, theta0, mu_E);
y0g = [a; e; i; OM; om; theta0];
y0 = [rr; vv];
th = thetaG(2024, 12, 22, 23, 20, 30,w_E);
groundTrack2(y0,t1,th,mu_E,w_E);
groundTrack2(y0,t2,th,mu_E,w_E);
groundTrack2(y0,t3,th,mu_E,w_E);

n = w_E*(k/m);
a_rep = power(mu_E/(n^2),1/3);
Tr = 2*pi*sqrt(a_rep^3/mu_E);
t1r = linspace(0, Tr, N);
t2r = linspace(0, 3*24*60*60, N);
[rrr, vvr] = kep2car(a_rep, e, i, OM, om, theta0, mu_E);
y0r = [rrr; vvr];
groundTrack2(y0r,t1r,th,mu_E,w_E);
groundTrack2(y0r,t2r,th,mu_E,w_E);
groundTrack2(y0r,t3,th,mu_E,w_E);

st = date2mjd2000(date);


fprintf('ne mancano solo 6\n');
fprintf('\n');

groundTrackJ2(y0,t1,th,mu_E,w_E, J2, R_E, mu_Moon, st);
groundTrackJ2(y0,t2,th,mu_E,w_E, J2, R_E, mu_Moon, st);
groundTrackJ2(y0,t3,th,mu_E,w_E, J2, R_E, mu_Moon, st);

fprintf('quasi arrivato\n');
fprintf('\n');

groundTrackJ2(y0r,t1r,th,mu_E,w_E, J2, R_E, mu_Moon, st);
groundTrackJ2(y0r,t2r,th,mu_E,w_E, J2, R_E, mu_Moon, st);
groundTrackJ2(y0r,t3,th,mu_E,w_E, J2, R_E, mu_Moon, st);

%% Le mucche
t3m = linspace(0, 30*T, N);
fprintf('Le mucche fanno mu, ma una fa mu Moon\n');
groundTrackJ2G(y0g,t1,th,mu_E,w_E, J2, R_E, mu_Moon, st);
groundTrackJ2G(y0g,t2,th,mu_E,w_E, J2, R_E, mu_Moon, st);
groundTrackJ2G(y0g,t3m,th,mu_E,w_E, J2, R_E, mu_Moon, st);