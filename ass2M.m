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
% LAST UPDATE: 30-12-2024
%
%--------------------------------------------------------------------------
% 2427 / Earth / J2 & Moon / 08-5-2057

%% SCRIPT INITIALISATION
clc
clear variables
close all

%% PRESENTATION
fprintf('----------------------------------------------\n');
fprintf('   PLANETARY EXPLORER MISSION ASSIGNMENT      \n');
fprintf('----------------------------------------------\n');
fprintf('\n');

%% IMPORTANT CONSTANTS
muE = astroConstants(13); % Earth's gravitational constant [km^3/s^2]
RE = astroConstants(23);  % Earth's radius [km]
wE = deg2rad(15.04)/3600; % Earth's angular velocity [rad/s]

muM = astroConstants(20); % Moon's gravitational constant [km^3/s^2]

J2 = astroConstants(9); % Earth's gravitational harmonic coefficient [-]

%% DATAS RECOVERY
% Below there are both the assigned values of the nominal orbit and the
% missing values (OM, om, theta and the starting date of propagation)
a = 23300; % Semi-major axis [km]
e = 0.523; % Eccentricity [-]
i = deg2rad(58.65); % Inclination [rad]
OM = 0; % Right ascension of ascending node [rad]
om = 0; % Pericenter anomaly [rad]
theta = 0; % True anomaly [rad]

date = [2057, 5, 8, 23, 20, 30]; % Starting date of propagation

k = 13; % Satellite's revolutions between ground track repetitions [-]
m = 3; % Earth's rotations between ground track repetition [-]

%% NOMINAL ORBIT
kep0 = [a, e, i, OM, om, theta]; % Keplerian parameters vector of nominal orbit

[r0, v0] = parorb2rv(a, e, i, OM, om, theta, muE); % Position and velocity vectors of initial point of nominal orbit
y0 = [r0; v0];

T = 2*pi*sqrt(a^3/muE); % Nominal orbit's period [s]

% Setting parameters for orbit representation
th0 = 0; % True anomaly of initial point [rad]
thf = 2*pi; % True anomaly of final point [rad]
dth = 0.1; % Resolution [rad]

Earth_3D
hold on
plot0rbit(a, e, i, OM, om, th0, thf, dth, muE, 2);
legend('-', 'Orbit', 'Pericenter', 'Eccentricity');
title('Nominal orbit');

%% GROUND TRACK OF NOMINAL NON PERTURBED ORBIT
% Computation of GST (thetaG0) for our simulation
jd = date2jd(date);
[J0, UT] = J0_computation(jd, date);
thetaG0 = thetaG0_computation(J0, UT, wE);

% Initialization of propagation's time intervals
N = 100000;
tv1 = linspace(0, T, N);
tv2 = linspace(0, 13*T, N);
tv3 = linspace(0, 30*T, N); % WE CAN CHANGE IT?

% Setting options for the ODE solver
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

% Propagation for the nominal non perturbed case
% T = 1
[~, Y1] = ode113 (@(t,y) ode_2bp(t,y,muE), tv1, y0, options);
Y1 = Y1';

[~, ~, ~, ~] = GroundTrack(Y1(1:3, :), thetaG0, tv1, wE);
title('Ground track for T = 1');

% T = 13
[~, Y2] = ode113 (@(t,y) ode_2bp(t,y,muE), tv2, y0, options);
Y2 = Y2';

[~, ~, ~, ~] = GroundTrack(Y2(1:3, :), thetaG0, tv2, wE);
title('Ground track for T = 13');

% T = 30
[~, Y3] = ode113 (@(t,y) ode_2bp(t,y,muE), tv3, y0, options); 
Y3 = Y3';

[~, ~, ~, ~] = GroundTrack(Y3(1:3, :), thetaG0, tv3, wE);
title('Ground track for T = 30');

%% REPEATING ORBIT CHARACTERISATION
% New semi-major axis to obtain repeating ground track
ar = a_repeatingGT(m, k, wE, muE); % [km]

[r0r, v0r] = parorb2rv (ar, e, i, OM, om, theta, muE); % Position and velocity vectors of initial point of repeating orbit
y0r = [r0r; v0r];

Tr =  2*pi*sqrt(ar^3/muE); % Repeating orbit's period [s]

Earth_3D
hold on
plot0rbit(ar, e, i, OM, om, th0, thf, dth, muE, 2);
legend('-', 'Orbit', 'Pericenter', 'Eccentricity');
title('Repeating orbit');

%% GROUND TRACK OF REPEATING NON PERTURBED ORBIT
% Initialization of propagation's time intervals
tv1r = linspace(0, Tr, N);
tv2r = linspace(0, 13*Tr, N);
tv3r = linspace(0, 30*Tr, N);

% Propagation for the repeating non perturbed case
% T = 1
[~, Y1r] = ode113 (@(t,y) ode_2bp(t,y,muE), tv1r, y0r, options);
Y1r = Y1r';

[~, ~, ~, ~] = GroundTrack(Y1r(1:3, :), thetaG0, tv1r, wE);
title('Repeating ground track for T = 1');

% T = 13
[~, Y2r] = ode113 (@(t,y) ode_2bp(t,y,muE), tv2r, y0r, options); 
Y2r = Y2r';

[~, ~, ~, ~] = GroundTrack(Y2r(1:3, :), thetaG0, tv2r, wE);
title('Repeating ground track for T = 13');

% T = 30
[~, Y3r] = ode113 (@(t,y) ode_2bp(t,y,muE), tv3r, y0r, options); 
Y3r = Y3r';

[~, ~, ~, ~] = GroundTrack(Y3r(1:3, :), thetaG0, tv3r, wE);
title('Repeating ground track for T = 30');

%% GROUND TRACK OF NOMINAL PERTURBED ORBIT
% Propagation for the nominal non perturbed case
% T = 1
[~, Y1p] = ode113 (@(t,y) ode_2bp_perturbed(t, y, muE, muM, RE, J2, date), tv1, y0, options); 
Y1p = Y1p'; 

[~, ~, ~, ~] = GroundTrack(Y1p(1:3, :), thetaG0, tv1, wE);
title('Ground track for perturbed orbit for T = 1');

% T = 13
[~, Y2p] = ode113 (@(t,y) ode_2bp_perturbed(t, y, muE, muM, RE, J2, date), tv2, y0, options);
Y2p = Y2p';

[~, ~, ~, ~] = GroundTrack(Y2p(1:3, :), thetaG0, tv2, wE);
title('Ground track for perturbed orbit for T = 13');

% T = 30
[~, Y3p] = ode113 (@(t,y) ode_2bp_perturbed(t, y, muE, muM, RE, J2, date), tv3, y0, options); 
Y3p = Y3p';

[~, ~, ~, ~] = GroundTrack(Y3p(1:3, :), thetaG0, tv3, wE);
title('Ground track for perturbed orbit for T = 30');