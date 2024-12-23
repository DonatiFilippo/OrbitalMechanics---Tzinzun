%% Main script for Interplanetary Mission assignment %%
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
% LAST UPDATE: 23-12-2024
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
[t_opt1, dv_min_ga] = ga(@(t) interplanetary1(t(1),t(2),t(3)), 3, [], [], [], [], lower, upper, [], options_ga);

% Usar el resultado del AG como punto inicial para fmincon
x0 = t_opt1; % t_opt es la solución inicial del AG
% Convert dates to Gregorian format
date_dep_ga = mjd20002date(t_opt1(1)); % Departure date
date_fb_ga = mjd20002date(t_opt1(2)); % Flyby date
date_arr_ga = mjd20002date(t_opt1(3)); % Arrival date


%% Refinement with FMINCON // Gradient based optimization 
% fmincon Configuration
options_fmincon = optimoptions('fmincon','Display', 'iter-detailed', 'Algorithm', 'sqp','StepTolerance', 1e-10, 'OptimalityTolerance', 1e-6);

fprintf('Refining Solution with FMINCON...\n');
[t_refined, dv_min_refined] = fmincon(@(t) interplanetary1(t(1), t(2), t(3)), t_opt1, [], [], [], [], lower, upper, [], options_fmincon);

% Convert refined dates to Gregorian format
date_dep_ref = mjd20002date(t_refined(1)); % Departure date
date_fb_ref = mjd20002date(t_refined(2)); % Flyby date
date_arr_ref = mjd20002date(t_refined(3)); % Arrival date

%% Algoritmo diferencial
options_de = optimoptions('particleswarm','SwarmSize', 200,'MaxIterations', 500,'Display', 'iter', 'UseParallel', true);

% Resolver con PSO / Particle Swarm
[t_refined2, dv_min_pso] = particleswarm(@(t) interplanetary1(t(1), t(2), t(3)),3, lower, upper, options_de);

date_dep_ref2 = mjd20002date(t_refined2(1)); % Departure date
date_fb_ref2 = mjd20002date(t_refined2(2)); % Flyby date
date_arr_ref2 = mjd20002date(t_refined2(3)); % Arrival date

%% Refinamiento Simulated Annealing

options_sa = optimoptions('simulannealbnd', 'MaxIterations', 2000,'Display', 'iter', 'PlotFcns', {@saplotbestx, @saplotbestf, @saplottemperature});

[t_refined_sa, dv_min_sa] = simulannealbnd(@(t) interplanetary1(t(1), t(2), t(3)), t_refined2, lower, upper, options_sa);

date_dep_sa = mjd20002date(t_refined_sa(1)); % Departure date
date_fb_sa = mjd20002date(t_refined_sa(2)); % Flyby date
date_arr_sa = mjd20002date(t_refined_sa(3)); % Arrival date

%% Resultados

fprintf('\n Genetic Algorithm Results:\n');
fprintf('Departure: %02d/%02d/%04d\n', date_dep_ga(3), date_dep_ga(2), date_dep_ga(1));
fprintf('Flyby: %02d/%02d/%04d\n', date_fb_ga(3), date_fb_ga(2), date_fb_ga(1));
fprintf('Arrival: %02d/%02d/%04d\n', date_arr_ga(3), date_arr_ga(2), date_arr_ga(1));
fprintf('Delta-v: %.2f km/s\n', dv_min_ga);


fprintf('\n FMINCON/Local Refinment Results:\n');
fprintf('Departure: %02d/%02d/%04d\n', date_dep_ref(3), date_dep_ref(2), date_dep_ref(1));
fprintf('Flyby: %02d/%02d/%04d\n', date_fb_ref(3), date_fb_ref(2), date_fb_ref(1));
fprintf('Arrival: %02d/%02d/%04d\n', date_arr_ref(3), date_arr_ref(2), date_arr_ref(1));
fprintf('Delta-v: %.2f km/s\n', dv_min_refined);


fprintf('\n PSO/Differential Algorithm Results:\n');
fprintf('Departure: %02d/%02d/%04d\n', date_arr_ref2(1));
fprintf('Flyby: %02d/%02d/%04d\n', date_arr_ref2(2));
fprintf('Arrival: %02d/%02d/%04d\n', date_arr_ref2(3));
fprintf('Delta-v: %02d/%02d/%04d\n', dv_min_pso);


fprintf('\n Simulated Annealing Results:\n');
fprintf('Departure: %02d/%02d/%04d \n', date_arr_sa(1));
fprintf('Flyby: %02d/%02d/%04d \n', date_arr_sa(2));
fprintf('Arrival: %02d/%02d/%04d \n', date_arr_sa(3));
fprintf('Delta-v: %02d/%02d/%04d \n', dv_min_sa);

%% GENERATE CHOP PLOT
fprintf('Generating CHOP Plot...\n');

% Convert MJD2000 to Gregorian dates for departure and arrival
t_dep_greg = arrayfun(@mjd20002date, t_dep, 'UniformOutput', false); % Convert each MJD2000 in t_dep
t_arr_greg = arrayfun(@mjd20002date, t_arr, 'UniformOutput', false); % Convert each MJD2000 in t_arr

% Convert the cell arrays to readable format for plotting
t_dep_labels = cellfun(@(x) datestr(datenum(x(1:3))), t_dep_greg, 'UniformOutput', false); % Convert to date strings
t_arr_labels = cellfun(@(x) datestr(datenum(x(1:3))), t_arr_greg, 'UniformOutput', false); % Convert to date strings


% Initialize delta-v matrix
delta_v = NaN(length(t_dep), length(t_arr));

% Calculate delta-v for all combinations of departure and arrival dates
for i = 1:length(t_dep)
    for j = 1:length(t_arr)
        if t_arr(j) > t_dep(i) % Ensure arrival is after departure
            % Call the interplanetary function
            delta_v(i, j) = interplanetary1(t_dep(i), (t_dep(i) + t_arr(j)) / 2, t_arr(j));
        end
    end
end


% Plot the CHOP Plot
figure;
hold on;
contourf(t_dep, t_arr, delta_v', 20, 'LineColor', 'none'); % Contour plot
colorbar; % Add a color bar
colormap('jet'); % Use jet colormap

% Customize the plot
xlabel('Departure Time');
ylabel('Arrival Time');
title('Delta-v Contour Plot');

% Update X and Y ticks with readable labels
set(gca, 'XTick', t_dep(1:30:end), 'XTickLabel', t_dep_labels(1:30:end));
set(gca, 'YTick', t_arr(1:30:end), 'YTickLabel', t_arr_labels(1:30:end));

% Highlight a specific chosen window
chosen_departure = t_dep(round(length(t_dep) / 2)); % Example chosen departure date
chosen_arrival = t_arr(round(length(t_arr) / 2));  % Example chosen arrival date
plot(chosen_departure, chosen_arrival, 'p', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'cyan'); % Add marker
legend('Chosen time window', 'Location', 'northeast');

hold off;