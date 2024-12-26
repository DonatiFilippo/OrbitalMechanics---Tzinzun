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

Re = astroConstants(23);

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
% We want to have an idea of the time it would take to complete the mission to 
% understand the ranges of available windows, that is, if the synodic period 
% between E-As is every 1.5 years, then the flight period between Mercury and 
% Earth must conclude before this period ends 
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

% Step size for iterating through time windows
step = 7; % We take a step of 1 day (adaptable)

% Calculate the synodic period with the most relevance
SP = max([T_syn_dep2fb, T_syn_fb2arr, T_syn_dep2arr]) / 86400; % Synodic period in days
SM = 0.4; % Safety margin for time of flight, 40% adjustment based on bibliography

% Time of flight ranges with safety margins
tof_t1_min = (1 - SM) * tof_t1; % Minimum time of flight Mercury -> Earth
tof_t1_max = (1 + SM) * tof_t1; % Maximum time of flight Mercury -> Earth

tof_t2_min = (1 - SM) * tof_t2; % Minimum time of flight Earth -> Asteroid
tof_t2_max = (1 + SM) * tof_t2; % Maximum time of flight Earth -> Asteroid

% Calculate the last possible departure time from Mercury
t_ldM = SP - tof_t1_min; % Last departure from Mercury to arrive within the synodic period

% Define the departure window from Mercury
w_dep = mdj_dep : step : mdj_dep + t_ldM; % First departure window from Mercury

% Calculate the arrival window at Earth
w_fb_min = w_dep(1) + tof_t1_min; % Earliest arrival at Earth
w_fb_max = w_dep(end) + tof_t1_max; % Latest arrival at Earth
w_fb = w_fb_min : step : w_fb_max; % Arrival window at Earth

% Calculate the departure window from Earth to the asteroid
w_arr_min = w_fb(1) + tof_t2_min; % Earliest departure from Earth to the asteroid
w_arr_max = w_fb(end) + tof_t2_max; % Ensure compatibility with the synodic period
w_arr = w_arr_min : step : w_arr_max; % Arrival window at the asteroid


%% Porkchop patches

days = 200;
dep_dates_1 = mdj_dep:days:(mdj_dep + t_ldM); % Fechas de salida para Mercury → Earth
arr_dates_1 = w_fb_min:days:w_fb_max;        % Fechas de llegada para Mercury → Earth
dep_dates_2 = w_arr_min:days:w_arr_max;      % Fechas de salida para Earth → Asteroid
arr_dates_2 = (mdj_arr - tof_t2_max):days:mdj_arr; % Fechas de llegada para Earth → Asteroid


% Cálculo de los deltas-v
[delta_v1, delta_v2] = porkchop(dep_dates_1, arr_dates_1, dep_dates_2, arr_dates_2);


[X1, Y1] = meshgrid(arr_dates_1, dep_dates_1); % Para delta_v1
[X2, Y2] = meshgrid(arr_dates_2, dep_dates_2); % Para delta_v2


figure()
hold on
grid on
title('Pork chop plot contour from Mercury to Earth')
xlabel('Time of arrival [MJD2000]')
ylabel('Time of departure [MJD2000]')
contour(X1, Y1, delta_v1, 50); % Usa la malla y matriz correspondientes
colorbar
colormap jet
datetick('x', 'dd/mm/yyyy', 'keepticks', 'keeplimits')
datetick('y', 'dd/mm/yyyy', 'keepticks', 'keeplimits')
set(gca, 'XTickLabelRotation', 45)
set(gca, 'YTickLabelRotation', 45)
hold off

figure()
hold on
grid on
title('Pork chop plot contour from Earth to Asteroid')
xlabel('Time of arrival [MJD2000]')
ylabel('Time of departure [MJD2000]')
contour(X2, Y2, delta_v2, 50);
colorbar
colormap jet
datetick('x', 'dd/mm/yyyy', 'keepticks', 'keeplimits')
datetick('y', 'dd/mm/yyyy', 'keepticks', 'keeplimits')
set(gca, 'XTickLabelRotation', 45)
set(gca, 'YTickLabelRotation', 45)
hold off




%% BEST SOLUTION FINDER ALGORITHMS
%% Genetic algorithm
lower = [w_dep(1) w_fb(1) w_arr(1)];           
upper = [w_dep(end) w_fb(end) w_arr(end)];               

lower_ga = [w_dep(1) w_fb(1) w_arr(1)];
upper_ga = [w_dep(end) w_fb(end) w_arr(end)];

% Options for genetic
options_ga = optimoptions('ga', 'PopulationSize', 500, ...
    'FunctionTolerance', 0.01, 'Display', 'off', 'MaxGenerations', 200);

% Solver
N = 1; % Number of departure windows examined
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
            lower_ga = lower;
            upper_ga = upper;
        end
        elapsedTime = toc(startTime);
        fprintf('Elapsed time : \n\n');
        fprintf('\n\n\n\n\n\n\n\n');
        fprintf('\b\b\b\b\b\b\b\b\b\b\b%6.2f s', elapsedTime);
        fprintf('\n\n');
    end

    lower = lower + t_ldM;
    upper = lower + t_ldM;
end

% Results with ga
date_dep_ga = mjd20002date(ceil(t_opt_ga(1)));
date_fb_ga = mjd20002date(ceil(t_opt_ga(2)));
date_arr_ga = mjd20002date(ceil(t_opt_ga(3)));

%% Refinement with FMINCON
% fmincon Configuration sqp selection options
options_fmincon = optimoptions('fmincon','Display', 'iter-detailed', 'Algorithm', 'sqp','StepTolerance', 1e-10, 'OptimalityTolerance', 1e-6);

% Fmincon solver
fprintf('Refining Solution with FMINCON...\n');
[t_refined, dv_min_refined] = fmincon(@(t) interplanetary(t(1), t(2), t(3)),  t_opt_ga, [], [], [], [], lower_ga, upper_ga, [], options_fmincon);

% Convert refined dates to Gregorian format
date_dep_ref = mjd20002date(t_refined(1));
date_fb_ref = mjd20002date(t_refined(2));
date_arr_ref = mjd20002date(t_refined(3));

%% Gradient refining method
% Options for gradient
options_grad = optimoptions('fminunc', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals', 1e4, 'MaxIter', 1e4, 'Display', 'off', 'Algorithm', 'quasi-newton'); 

% Gradient solver
[t_opt_grad, dv_min_grad] = fminunc(@(t) interplanetary(t(1), t(2), t(3)), t_opt_ga, options_grad);

% Results with gradient
date_dep_grad = mjd20002date(ceil(t_opt_grad(1)));
date_fb_grad = mjd20002date(ceil(t_opt_grad(2)));
date_arr_grad = mjd20002date(ceil(t_opt_grad(3)));

%% Refinamiento Simulated Annealing

options_sa = optimoptions('simulannealbnd', 'MaxIterations', 2000,'Display', 'iter', 'PlotFcns', {@saplotbestx, @saplotbestf, @saplottemperature});

[t_refined_sa, dv_min_sa] = simulannealbnd(@(t) interplanetary(t(1), t(2), t(3)), t_opt_ga, lower_ga, upper_ga, options_sa);

date_dep_sa = mjd20002date(t_refined_sa(1));
date_fb_sa = mjd20002date(t_refined_sa(2));
date_arr_sa = mjd20002date(t_refined_sa(3));

%% Algorithm comparison
% Genetic Algorithm
fprintf('\n\n')
fprintf('Genetic Algorithm Results:\n\n');
fprintf('Departure: %02d/%02d/%04d\n', date_dep_ga(3), date_dep_ga(2), date_dep_ga(1));
fprintf('Flyby: %02d/%02d/%04d\n', date_fb_ga(3), date_fb_ga(2), date_fb_ga(1));
fprintf('Arrival: %02d/%02d/%04d\n', date_arr_ga(3), date_arr_ga(2), date_arr_ga(1));
fprintf('Delta-v: %.2f km/s\n', dv_min_ga);
fprintf('\n\n')

% FMINCON Refinement
fprintf('FMINCON/Local Refinement Results:\n\n');
fprintf('Departure: %02d/%02d/%04d\n', date_dep_ref(3), date_dep_ref(2), date_dep_ref(1));
fprintf('Flyby: %02d/%02d/%04d\n', date_fb_ref(3), date_fb_ref(2), date_fb_ref(1));
fprintf('Arrival: %02d/%02d/%04d\n', date_arr_ref(3), date_arr_ref(2), date_arr_ref(1));
fprintf('Delta-v: %.2f km/s\n', dv_min_refined);
fprintf('\n\n')

% Refined solution with gradient
fprintf('Refined solution with gradient :\n\n');
fprintf('Departure date : %02d/%02d/%04d\n', date_dep_grad(3), date_dep_grad(2), date_dep_grad(1));
fprintf('Fly-by date: %02d/%02d/%04d\n', date_fb_grad(3), date_fb_grad(2), date_fb_grad(1));
fprintf('Arrival date : %02d/%02d/%04d\n', date_arr_grad(3), date_arr_grad(2), date_arr_grad(1));
fprintf('Minimised cost with gradient : %f km/s \n', dv_min_grad);
fprintf('\n\n')

% Simulated Annealing
fprintf('Simulated Annealing Results:\n\n');
fprintf('Departure: %02d/%02d/%04d \n', date_dep_sa(3), date_dep_sa(2), date_dep_sa(1));
fprintf('Flyby: %02d/%02d/%04d \n', date_fb_sa(3), date_fb_sa(2), date_fb_sa(1));
fprintf('Arrival: %02d/%02d/%04d \n', date_arr_sa(3), date_arr_sa(2), date_arr_sa(1));
fprintf('Delta-v: %.2f km/s\n', dv_min_sa);
fprintf('\n\n')

