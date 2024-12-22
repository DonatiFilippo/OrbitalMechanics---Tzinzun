%% Main script for Interplanetary Mission assignment %%
%
% Group number : 27
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
% LAST UPDATE: 21-12-2024
%

% 2427 / Mercury / Earth / Asteroid N.30 / 00:00:00 01/01/2030 / 00:00:00 01/01/2060

%% SCRIPT INITIALISATION
clc
clear variables
close all

%% PRESENTATION
fprintf('----------------------------------------------\n');
fprintf('     INTERPLANETARY MISSION ASSIGNMENT        \n');
fprintf('----------------------------------------------\n');
fprintf('\n');

%% IMPOSED DATES
Departure_date = [2030, 1, 1, 0, 0, 0]; % Earliest departure date [Gregorian date]
Arrival_date = [2060, 1, 1, 0, 0, 0]; % Latest arrival date [Gregorian date]
 
mdj_dep = date2mjd2000(Departure_date); % Earliest departure date conversion from Gregorian date to modified Julian day 2000 number
mdj_arr = date2mjd2000(Arrival_date); % Latest arrival date conversion from Gregorian date to modified Julian day 2000 number

%% PHYSICAL PARAMETERS
Departure_planet = 1; % Mercury as the departure planet 
Flyby_planet = 3; % Earth as the flyby planet
Arrival_asteroid_id = 30; % Asteroid no.30 as the arrival objective

muM = astroConstants(11); % Mercury's gravitational constant [km^3/s^2]
muE = astroConstants(13); % Earth's gravitational constant [km^3/s^2]
muS = astroConstants(4);  % Sun's gravitational constant [km^3/s^2]

[kep_dep, ~] = uplanet(0, Departure_planet); % Mercury's keplerian elements at initial time / kep = [a e i Om om theta] [km, rad]
[kep_fb, ~] = uplanet(0, Flyby_planet); % Earth's keplerian elements at initial time / kep = [a e i Om om theta] [km, rad]
[kep_arr, ~, ~] = ephAsteroids(0, Arrival_asteroid_id); % Asteroid's keplerian elements at initial time / kep = [a e i Om om theta] [km, rad]

%% PERIODS COMPUTING
T_dep = 2*pi*sqrt(kep_dep(1)^3/muS); % Mercury's orbital period [s]
T_fb = 2*pi*sqrt(kep_fb(1)^3/muS); % Earth's orbital period [s]
T_arr = 2*pi*sqrt(kep_arr(1)^3/muS); % Asteroid's orbital period [s]

T_syn_dep2fb = T_dep * T_fb/abs(T_dep - T_fb);  % Mercury's synodic orbital period with respect to Earth [s]
T_syn_fb2arr = T_fb * T_arr/abs(T_fb - T_arr); % Earth's synodic orbital period with respect to the asteroid [s]
T_syn_dep2arr = T_dep * T_arr/abs(T_dep - T_arr); % Mercury's synodic orbital period with respect to the asteroid [s]

fprintf('Displaying synodics period in years :\n\n');
fprintf('Mercury''s synodic orbital period with respect to Earth : %f years \n', T_syn_dep2fb/(86400*365.25));
fprintf('Earth''s synodic orbital period with respect to the asteroid : %f years \n', T_syn_fb2arr/(86400*365.25));
fprintf('Mercury''s synodic orbital period with respect to the asteroid : %f years \n', T_syn_dep2arr/(86400*365.25));
fprintf('\n\n');

% The greatest synodic period is the one of Earth with respect to the asteroid
% To find the best window we have to search with respect to this period
% We will admit that this synodic period is 1.5 years long to simplify

%% WINDOWS APPROXIMATION
% We want to have an idea of the time it would take to do the mission
% We can compute the time of flight for Hohmann transfers to get a first view
% Orbits will be considered coplanar and circular

a_t1 = (kep_dep(1) + kep_fb(1))/2; % Semi-major axis of the first transfert arc between Mercury and Earth [km]
a_t2 = (kep_fb(1) + kep_arr(1))/2; % Semi-major axis of the second transfert arc between Earth and the asteroid [km]

T_t1 = 2*pi*sqrt(a_t1^3/muS); % Period of the first transfer arc [s]
T_t2 = 2*pi*sqrt(a_t2^3/muS); % Period of the second transfer arc [s]

tof_t1 = 1/2*T_t1 / 86400; % Time of flight of the first transfer arc [days]
tof_t2 = 1/2*T_t2 / 86400; % Time of flight of the second transfer arc [days]
tof_h = (tof_t1 + tof_t2); % Total time of flight for Hohmann transfer [days]

fprintf('Displaying tof of transfer arc for Hohmann transfers in years :\n\n');
fprintf('Time of flight of the first transfer arc from Mercury to Earth : %f years\n', tof_t1/(365.25));
fprintf('Time of flight of the second transfer arc from Earth to the asteroid : %f years \n', tof_t2/(365.25));
fprintf('Total time of flight from Mercury to the asteroid : %f years \n', tof_h/(365.25));
fprintf('\n\n');

% The real time of flight cannot be computed for now
% Thus, we have to rely on the Hohmann's one to approximate the real one
% Recall : the synodic period of Earth wrt the asteroid is 1.5 years
% The Hohamnn tof to go from Mercury to Earth is approximately 4 months
% The Hohamnn tof to go from Earth the asteroid is approximately 13 months
% We consider the scenarios where it takes 50% to 150% of these tof 

dt_dep = 1.5*365.25; % Duration of the first window with respect to the synodic periods (2 years) [days]
step = 1; % We take a step of 1 day (adaptable)

t_dep = mdj_dep : step : mdj_dep + dt_dep; % First departure window from Mercury
t_fb = t_dep(1) + 1/2*tof_t1 : step : t_dep(end) + 3/2*tof_t1; % First arrival window on Earth
t_arr = t_fb(1) + 1/2*tof_t2 : step : t_fb(end) + 3/2*tof_t2; % First arrival window on the asteroid

% We now consider that we should find an optimal solution between t_dep(1) and t_arr(end)

%% BEST SOLUTION FINDER ALGORITHMS
% Genetic algorithm
lower = [t_dep(1) t_fb(1) t_arr(1)];           
upper = [t_dep(end) t_fb(end) t_arr(end)];               

% Options for genetic
options_ga = optimoptions('ga', 'PopulationSize', 100, 'Display', 'iter', 'MaxGenerations', 100);

% Solving
[t_opt, dv_min] = ga(@(t) interplanetary(t(1),t(2),t(3)), 3, [], [], [], [], lower, upper, [], options_ga);

% Print results
fprintf('Résultats optimisés avec ga :\n');
fprintf('Départ : %.2f jours julien\n', t_opt(1));
fprintf('Survol : %.2f jours julien\n', t_opt(2));
fprintf('Arrivée : %.2f jours julien\n', t_opt(3));
fprintf('Delta-v total : %.2f km/s\n', dv_min);