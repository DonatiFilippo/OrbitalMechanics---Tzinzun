%% Main script for Planetary Explorer Mission assignment %% LOCAL VERSION
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

% 2427 / Earth / J2 & Moon / 22-12-2024  

%% SCRIPT INITIALISATION
clc
clear variables
close all

%% PRESENTATION
fprintf('----------------------------------------------\n');
fprintf('   PLANETARY EXPLORER MISSION ASSIGNMENT      \n');
fprintf('----------------------------------------------\n');
fprintf('\n');

%% GENERAL DATAS
muE = astroConstants(13);
RE = astroConstants(23);
wE = deg2rad(15.04)/3600;

muM = astroConstants(20);

J2 = astroConstants(9);

%% ORBIT CHARACTERISATION
a = 23300; % Semi-major axis [km]
e = 0.523; % Eccentricity [-]
i = deg2rad(58.65); % Inclination [rad]
OM = 0; % RAAN [rad]
om = 0;
theta = 0;
N = 100000;

kep0 = [a, e, i, OM, om, theta];
[r0, v0] = parorb2rv(a, e, i, OM, om, theta, muE);
%%
y0 = [r0; v0];
T = 2*pi*sqrt(a^3/muE); % Orbital period [s]

% Orbit representation
% Earth_3D
% hold on
% plot0rbit(a, e, i, OM, om, 0, 2*pi, 0.1, muE, 2)

% ThetaG0 computation
date = [2024, 12, 22, 23, 20, 30];
jd = date2jd(date);
[J0, UT] = J0_computation(jd, date);
thetaG0 = thetaG0_computation(J0, UT, wE);

%% Groud track not repeating
t0 = 0;
tv = linspace(t0, 50*T, N);

% Options for the ODE solver
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

% Perform the integration for the standard (non perturbed) case
[Tout, Yout] = ode113 (@(t,y) ode_2bp(t,y,muE), tv, y0, options);
Tout = Tout'; % [1xl] row vector
Yout = Yout'; % [6xl] row vector

[Toutp, Youtp] = ode113 (@(t,y) ode_2bp_perturbed(t, y, muE, muM, RE, J2, date), tv, y0, options);
Toutp = Toutp'; % [1xl] row vector
Youtp = Youtp'; % [6xl] row vector

for k = 1:length(tv)
    [a_car(k), e_car(k), i_car(k), OM_car(k), om_car(k), theta_car(k)] = rv2parorb(Youtp(1:3, k), Youtp(4:end, k), muE);
end


% aspetta1 = GroundTrack(Yout(1:3, :), thetaG0, tv, wE);
% 
% dontneed = GroundTrack(Youtp(1:3, :), thetaG0, tv, wE);

%% Orbit propagation using Gauss
[ToutKEP, KEPout] = ode113 (@(t,kep)  Gauss_TNH(t, kep, muE, RE, J2, muM, date), tv, kep0, options);
KEPout = KEPout';

% [ToutKEP_RSW, KEPout_RSW] = ode113 (@(t,kep)  Gauss_RSW(t, kep, muE, RE, J2, muM, date), tv, kep0, options);
% KEPout_RSW = KEPout_RSW';
% 
% whoareu = GroundTrackKEP(KEPout, thetaG0, tv, wE, muE);

%% %% Repeating ground track
% k = 13;
% m = 3;
% ar = a_repeatingGT(m, k, wE, muE);
% 
% Tr =  2*pi*sqrt(ar^3/muE);
% tvr = linspace(t0, 10*24*60*60, 10000);
% 
% [r0r, v0r] = parorb2rv (ar, e, i, OM, om, theta0, muE);
% y0r = [r0r; v0r];
% 
% [Toutr, Youtr] = ode113 (@(t,y) ode_2bp(t,y,muE), tvr, y0r, options);
% Toutr = Toutr'; % [1xl] row vector
% Youtr = Youtr'; % [6xl] row vector
% 
% [Tout_rp, Yout_rp] = ode113 (@(t,y) ode_2bp_perturbed(t, y, muE, muM, RE, J2, date), tvr, y0r, options);
% Tout_rp = Tout_rp'; % [1xl] row vector
% Yout_rp = Yout_rp'; % [6xl] row vector
% 
% [alpha, delta, lon, lat] = GroundTrack(Youtr(1:3, :), thetaG0, tvr, wE);
% [alpha, delta, lon, lat] = GroundTrack(Yout_rp(1:3, :), thetaG0, tvr, wE);

%% PART 4: PLOT OF KEPLERIAN ELEMENTS
figure
plot(ToutKEP'/T, KEPout(1, :), 'g-')
grid on
title("semi-major axis")

figure
semilogy(Tout/T, abs(a_car-KEPout(1,:)), 'b-')
grid on
title("semi-major axis error")

figure
plot(ToutKEP'/T, KEPout(2, :), 'g-')
grid on
title("eccentricity")

figure
semilogy(Tout/T, abs(e_car-KEPout(2,:)), 'b-')
grid on
title("eccentricity error")

figure
plot(ToutKEP'/T, KEPout(3, :), 'g-')
grid on
title("inclination")

figure
semilogy(Tout/T, abs(i_car-KEPout(3,:)), 'b-')
grid on
title("inclination error")

figure
plot(ToutKEP'/T, KEPout(4, :), 'g-')
grid on
title("RAAN")

figure
semilogy(Tout/T, abs(OM_car-KEPout(4,:)), 'b-')
grid on
title("RAAN error")

figure
plot(ToutKEP'/T, KEPout(5, :), 'g-')
grid on
title("perigee anomaly")

figure
semilogy(Tout/T, abs(om_car-KEPout(5,:)), 'b-')
grid on
title("perigee anomaly error")

figure
plot(ToutKEP'/T, KEPout(6, :), 'g-')
grid on
title("true anomaly anomaly")

figure
semilogy(Tout/T, abs(unwrap(theta_car)-KEPout(6,:)), 'b-')
grid on
title("true anomaly error error")

%% 5 PLOT PERTURBED ORBIT

%% 6 FILTERING (questa parte devo farla meglio perchè avevo cambaito una cosa ed esce storta)
N = length(tv);

a_lt = movmean(KEPout(1,:), N/500);
figure
plot(ToutKEP'/T, KEPout(1,:), 'g-', ToutKEP'/T, a_lt, 'b-')
grid on


e_lt = movmean(KEPout(2,:), N/500);
e_s = movmean(KEPout(2,:), N/15.4010);
figure
plot(ToutKEP'/T, KEPout(2,:), 'g-', ToutKEP'/T, e_lt, 'b-', ToutKEP'/T, e_s, 'r-')
grid on


i_lt = movmean(KEPout(3,:), N/500);
i_s = movmean(KEPout(3,:), N/15.1699);
figure
plot(ToutKEP'/T, KEPout(3,:), 'g-', ToutKEP'/T, i_lt, 'b-', ToutKEP'/T, i_s, 'r-')
grid on


OM_lt = movmean(KEPout(4,:), N/50);
%OM_s = movmean(KEPout(4,:), N/15.1699);
figure
plot(ToutKEP'/T, KEPout(4,:), 'g-', ToutKEP'/T, OM_lt, 'b-')
grid on
