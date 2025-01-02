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
legend('', 'Orbit', 'Pericenter', 'Eccentricity');
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
N = 100000;
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
% Propagation for the nominal perturbed case
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

%% GROUND TRACK OF REPEATING PERTURBED ORBIT
% Propagation for the repeating perturbed case
% T = 1
[~, Y1rp] = ode113 (@(t,y) ode_2bp_perturbed(t, y, muE, muM, RE, J2, date), tv1r, y0r, options); 
Y1rp = Y1rp'; 

[~, ~, ~, ~] = GroundTrack(Y1rp(1:3, :), thetaG0, tv1r, wE);
title('Repeating ground track for perturbed orbit for T = 1');

% T = 13
[~, Y2rp] = ode113 (@(t,y) ode_2bp_perturbed(t, y, muE, muM, RE, J2, date), tv2r, y0r, options);
Y2rp = Y2rp';

[~, ~, ~, ~] = GroundTrack(Y2rp(1:3, :), thetaG0, tv2r, wE);
title('Repeating ground track for perturbed orbit for T = 13');

% T = 30
[~, Y3rp] = ode113 (@(t,y) ode_2bp_perturbed(t, y, muE, muM, RE, J2, date), tv3r, y0r, options); 
Y3rp = Y3rp';

[~, ~, ~, ~] = GroundTrack(Y3rp(1:3, :), thetaG0, tv3r, wE);
title('Repeating ground track for perturbed orbit for T = 30');

%% PROPAGATION OF NOMINAL PERTURBED ORBIT WITH DIFFERENT METHODS
% Setting the time span for the propagation
N = 100000;
tv = linspace(0, 100*T, N);

% Propagation in Cartesian coordinates
[TC, Y] = ode113 (@(t,y) ode_2bp_perturbed(t, y, muE, muM, RE, J2, date), tv, y0, options); 
Y = Y';
TC = TC';

% Keplerian elements recovery
for k = 1:N
    [a_car(k), e_car(k), i_car(k), OM_car(k), om_car(k), theta_car(k)] = rv2parorb(Y(1:3, k), Y(4:end, k), muE);
end
OM_car = unwrap(OM_car);
om_car = unwrap(om_car);
theta_car = unwrap(theta_car);

% Propagation of Keplerian elements using Gauss' equations
[TKEP, KEP] = ode113 (@(t,kep)  Gauss_TNH(t, kep, muE, RE, J2, muM, date), tv, kep0, options);
KEP = KEP';
TKEP = TKEP';

%% PROPAGATION METHODS COMPARISON
% Conversion to degrees of Keplerian elements obtained from Cartesian
% coordinates propagation for readability
i_car = rad2deg(i_car);
OM_car = rad2deg(OM_car);
om_car = rad2deg(om_car);
theta_car = rad2deg(theta_car);

% Conversion to degrees of Keplerian elements obtained from Gauss'
% equations propagation
KEP(3:6, :) = rad2deg(KEP(3:6, :));

% Conversion to degrees of nominal Keplerian elements
kep0(3:6) = rad2deg(kep0(3:6));

% Evolutions over time of Keplerian elements obtained with both methods and
% comparison
% -------------------------------------------------------------------------

% Semi-major axis
figure
plot(TC/T, a_car, 'b-')
grid on
xlabel('Time [T]')
ylabel('a_c_a_r [km]')
title("Evolution of semi-major axis Cartesian method")

figure
plot(TKEP/T, KEP(1,:), 'b-')
grid on
xlabel('Time [T]')
ylabel('a_g_a_u_s_s [km]')
title("Evolution of semi-major axis Gauss method")

figure
plot(TKEP/T, abs(a_car-KEP(1,:))/kep0(1), 'b-')
grid on
xlabel('Time [T]')
ylabel('|a_c_a_r - a_g_a_u_s_s|/a0 [-]')
title("Semi-major axis' relative error")

% -------------------------------------------------------------------------

% Eccentricity
figure
plot(TC/T, e_car, 'b-')
grid on
xlabel('Time [T]')
ylabel('e_c_a_r [-]')
title("Evolution of eccentricity Cartesian method")

figure
plot(TKEP/T, KEP(2,:), 'b-')
grid on
xlabel('Time [T]')
ylabel('e_g_a_u_s_s [-]')
title("Evolution of eccentricity Gauss method")

figure
plot(TKEP/T, abs(e_car-KEP(2,:)), 'b-')
grid on
xlabel('Time [T]')
ylabel('|e_c_a_r - e_g_a_u_s_s| [-]')
title("Eccentricity's absolute error")

% -------------------------------------------------------------------------

% Inclination
figure
plot(TC/T, i_car, 'b-')
grid on
xlabel('Time [T]')
ylabel('i_c_a_r [deg]')
title("Evolution of inclination Cartesian method")

figure
plot(TKEP/T, KEP(3,:), 'b-')
grid on
xlabel('Time [T]')
ylabel('i_g_a_u_s_s [deg]')
title("Evolution of inclination Gauss method")

figure
plot(TKEP/T, abs(i_car-KEP(3,:))/(2*pi), 'b-')
grid on
xlabel('Time [T]')
ylabel('|i_c_a_r - i_g_a_u_s_s|/2π [-]')
title("Inclination's relative error")

% -------------------------------------------------------------------------

% Right ascension of ascending node
figure
plot(TC/T, OM_car, 'b-')
grid on
xlabel('Time [T]')
ylabel('Ω_c_a_r [deg]')
title("Evolution of RAAN Cartesian method")

figure
plot(TKEP/T, KEP(4,:), 'b-')
grid on
xlabel('Time [T]')
ylabel('Ω_g_a_u_s_s [deg]')
title("Evolution of RAAN Gauss method")

figure
plot(TKEP/T, abs(OM_car-KEP(4,:))/(2*pi), 'b-')
grid on
xlabel('Time [T]')
ylabel('|Ω_c_a_r - Ω_g_a_u_s_s|/2π [-]')
title("RAAN's relative error")

% -------------------------------------------------------------------------

% Pericenter anomaly
figure
plot(TC/T, om_car, 'b-')
grid on
xlabel('Time [T]')
ylabel('ω_c_a_r [deg]')
title("Evolution of pericenter anomaly Cartesian method")

figure
plot(TKEP/T, KEP(5,:), 'b-')
grid on
xlabel('Time [T]')
ylabel('ω_g_a_u_s_s [deg]')
title("Evolution of pericenter anomaly Gauss method")

figure
plot(TKEP/T, abs(om_car-KEP(5,:))/(2*pi), 'b-')
grid on
xlabel('Time [T]')
ylabel('|ω_c_a_r - ω_g_a_u_s_s|/2π [-]')
title("Pericenter anomaly's relative error")

% -------------------------------------------------------------------------

% True anomaly
figure
plot(TC/T, theta_car, 'b-')
grid on
xlabel('Time [T]')
ylabel('θ_c_a_r [deg]')
title("Evolution of true anomaly Cartesian method")

figure
plot(TKEP/T, KEP(6,:), 'b-')
grid on
xlabel('Time [T]')
ylabel('θ_g_a_u_s_s [deg]')
title("Evolution of true anomaly Gauss method")

figure
plot(TKEP/T, abs(theta_car-KEP(6,:))/(2*pi), 'b-')
grid on
xlabel('Time [T]')
ylabel('|θ_c_a_r - θ_g_a_u_s_s|/2π [-]')
title("Pericenter anomaly's relative error")

%% ORBIT'S EVOLUTION
% ATTENZIONE: LA RAPPRESENTAZIONE DELL'ORBITA PUO' ESSERE MIGLIORATA
% Setting the time span for the propagation
n_orbits = 1000;
tv_plot = linspace(0, n_orbits*T, N);

% Propagation in Cartesian coordinates
[T_plot, Y_plot] = ode113 (@(t,y) ode_2bp_perturbed(t, y, muE, muM, RE, J2, date), tv_plot, y0, options); 
Y_plot = Y_plot';
T_plot = T_plot';


Earth_3D
hold on
scatter3(Y_plot(1,:), Y_plot(2,:), Y_plot(3,:), 5, T_plot./T, 'filled');
colormap;
a = colorbar;
a.Title.String = "Orbital period [T]";
clim([1 n_orbits]);
xlabel('X [km]'); 
ylabel('Y [km]');
zlabel('Z [km]');
title('Perturbed two-body problem orbit');
axis equal;
grid on;

%% FILTERING
% 