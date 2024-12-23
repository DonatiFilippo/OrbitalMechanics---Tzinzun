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
%--------------------------------------------------------------------------

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
 
mdj_dep = date2mjd2000(Departure_date); % Earliest departure date converted to modified Julian day 2000
mdj_arr = date2mjd2000(Arrival_date); % Latest arrival date converted to modified Julian day 2000

%% PHYSICAL PARAMETERS
Departure_planet = 1; % Mercury as the departure planet 
Flyby_planet = 3; % Earth as the flyby planet
Arrival_asteroid_id = 30; % Asteroid no.30 as the arrival objective

muM = astroConstants(11); % Mercury's gravitational constant [km^3/s^2]
muE = astroConstants(13); % Earth's gravitational constant [km^3/s^2]
muS = astroConstants(4); % Sun's gravitational constant [km^3/s^2]

% Keplerian elements computing / kep = [a e i Om om theta] [km, rad]
[kep_dep, ~] = uplanet(0, Departure_planet); % Mercury's keplerian elements at initial time
[kep_fb, ~] = uplanet(0, Flyby_planet); % Earth's keplerian elements at initial time
[kep_arr, ~, ~] = ephAsteroids(0, Arrival_asteroid_id); % Asteroid's keplerian elements at initial time

%% PERIODS COMPUTING
T_dep = 2*pi*sqrt(kep_dep(1)^3/muS); % Mercury's orbital period [s]
T_fb = 2*pi*sqrt(kep_fb(1)^3/muS); % Earth's orbital period [s]
T_arr = 2*pi*sqrt(kep_arr(1)^3/muS); % Asteroid's orbital period [s]

T_syn_dep2fb = T_dep * T_fb/abs(T_dep - T_fb); % Mercury's synodic orbital period with respect to Earth [s]
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

%% WINDOWS RESEARCH
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

dt_windows = 1.5*365.25; % Duration of the first window with respect to the synodic periods (1.5 years) [days]
step = 1; % We take a step of 1 day (adaptable)

t_dep = mdj_dep : step : mdj_dep + dt_windows; % First departure window from Mercury
t_fb = t_dep(1) + 1/2*tof_t1 : step : t_dep(end) + 3/2*tof_t1; % First arrival window on Earth
t_arr = t_fb(1) + 1/2*tof_t2 : step : t_fb(end) + 3/2*tof_t2; % First arrival window on the asteroid

% We now consider that we should find an optimal solution between t_dep(1) and t_arr(end)

%% BEST SOLUTION FINDER ALGORITHMS
%% Genetic algorithm
lower = [t_dep(1) t_fb(1) t_arr(1)];           
upper = [t_dep(end) t_fb(end) t_arr(end)];               

% Options for genetic
options_ga = optimoptions('ga', 'PopulationSize', 500, ...
    'FunctionTolerance', 0.01, 'Display', 'off', 'MaxGenerations', 200);

% Solver
N = 18; % Number of departure windows examined
N_ga = 5; % Number of genetic algorithm iteration to have better results
dv_min_ga = 50; % Arbitrary chosen value of total cost
t_opt_ga = [0, 0, 0]; % Storage value for the chosen windows

fprintf('Genetic algorithm computing ... \n\n');
startTime = tic;
for i = 1:N
    fprintf('ITERATION NUMBER : %2.f \n \n', i);

    for j = 1:N_ga
        [t_opt_ga_computed, dv_min_ga_computed] = ga(@(t) interplanetary(t(1),t(2),t(3)), 3, [], [], [], [], lower, upper, [], options_ga);

        if dv_min_ga_computed < dv_min_ga
            dv_min_ga = dv_min_ga_computed;
            t_opt_ga = t_opt_ga_computed;
        end

        elapsedTime = toc(startTime);
        fprintf('Elapsed time : \n\n');
        fprintf('\n\n\n\n\n\n\n\n');
        fprintf('\b\b\b\b\b\b\b\b\b\b\b%6.2f s', elapsedTime);
        fprintf('\n\n');
    end

    lower = lower + dt_windows;
    upper = lower + dt_windows;
end
 
% Results with ga
dep_grad = mjd20002date(ceil(t_opt_ga(1)));
fb_grad = mjd20002date(ceil(t_opt_ga(2)));
arr_grad = mjd20002date(ceil(t_opt_ga(3)));

dep_date = sprintf('%02d/%02d/%04d %02d:%02d:%02d', ...
    dep_grad(2), dep_grad(3), dep_grad(1), dep_grad(4), dep_grad(5), dep_grad(6));
fb_date = sprintf('%02d/%02d/%04d %02d:%02d:%02d', ...
    fb_grad(2), fb_grad(3), fb_grad(1), fb_grad(4), fb_grad(5), fb_grad(6));
arr_date = sprintf('%02d/%02d/%04d %02d:%02d:%02d', ...
    arr_grad(2), arr_grad(3), arr_grad(1), arr_grad(4), arr_grad(5), arr_grad(6));

% Print  results with ga
fprintf('Optmised results with ga :\n\n');
fprintf('Departure date : %s \n', dep_date);
fprintf('Fly-by date: %s \n', fb_date);
fprintf('Arrival date : %s \n', arr_date);
fprintf('Total optimised cost (dv) : %f km/s \n', dv_min_ga);
fprintf('\n\n');

%% Gradient refining method
% Options for gradient
options_grad = optimoptions('fminunc', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals', 1e4, 'MaxIter', 1e4, 'Display', 'off', 'Algorithm', 'quasi-newton'); 

% Solver
[t_opt_grad, dv_min_grad] = fminunc(@(t) interplanetary(t(1), t(2), t(3)), t_opt_ga, options_grad);

% Results with gradient
dep_grad = mjd20002date(ceil(t_opt_grad(1)));
fb_grad = mjd20002date(ceil(t_opt_grad(2)));
arr_grad = mjd20002date(ceil(t_opt_grad(3)));

dep_date = sprintf('%02d/%02d/%04d %02d:%02d:%02d', ...
    dep_grad(2), dep_grad(3), dep_grad(1), dep_grad(4), dep_grad(5), dep_grad(6));
fb_date = sprintf('%02d/%02d/%04d %02d:%02d:%02d', ...
    fb_grad(2), fb_grad(3), fb_grad(1), fb_grad(4), fb_grad(5), fb_grad(6));
arr_date = sprintf('%02d/%02d/%04d %02d:%02d:%02d', ...
    arr_grad(2), arr_grad(3), arr_grad(1), arr_grad(4), arr_grad(5), arr_grad(6));

% Print refined solution with gradient
fprintf('Refined solution with gradient :\n\n');
fprintf('Departure date : %s \n', dep_date);
fprintf('Fly-by date: %s \n', fb_date);
fprintf('Arrival date : %s \n', arr_date);
fprintf('Minimised cost with gradient : %f \n', dv_min_grad);

