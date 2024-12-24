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
fprintf('Mercury''s synodic orbital period with respect to Earth : %f dias \n', T_syn_dep2fb/(86400*365.25)*365.25);
fprintf('Earth''s synodic orbital period with respect to the asteroid : %f dias \n', T_syn_fb2arr/(86400*365.25)*365.25);
fprintf('Mercury''s synodic orbital period with respect to the asteroid : %f days \n', T_syn_dep2arr/(86400*365.25)*365.25);
fprintf('\n\n');

% The greatest synodic period is the one of Earth with respect to the
% asteroid, se debe considerar princpipalmente este periodo porque es el
% que más tarda en volver a ser posible
% To find the best window we have to search with respect to this period
% We will admit that this synodic period is 1.5 years long to simplify

%% WINDOWS APPROXIMATION
% We want to have an idea of the time it would take to do the mission para
% conocer los rangos de ventanas disponibles, es decir, si el synodic
% period de E-As es cada 1.5 años, entonces el periodo de vuelo entre Mer y
% Ea tiene que concluir antes de que este periodo termine
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
fprintf('Time of flight of the first transfer arc from Mercury to Earth : %f days\n', tof_t1);
fprintf('Time of flight of the second transfer arc from Earth to the asteroid : %f days \n', tof_t2);
fprintf('Total time of flight from Mercury to the asteroid : %f days \n', tof_h);
fprintf('\n\n');

% The real time of flight cannot be computed for now
% Thus, we have to rely on the Hohmann's one to approximate the real one
% Recall : the synodic period of Earth wrt the asteroid is 1.5 years
% The Hohamnn tof to go from Mercury to Earth is approximately 4 months,
% así que la ventana de lanzamiento desde Mer debe terminar antes de este
% periodo para llegar a tiempo al synodic period
% The Hohamnn tof to go from Earth the asteroid is approximately 13 months
% We consider the scenarios where it takes 50% to 150% of these tof 


% Step size for iterating through time windows
step = 1; % We take a step of 1 day (adaptable)

% Calculate the synodic period with the most relevance
SP = max([T_syn_dep2fb, T_syn_fb2arr, T_syn_dep2arr]) / 86400; % Synodic period in days
SM = 0.4; % Safety margin for time of flight, 40% adjustment based on bibliography I put on Latex file


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


%% BEST SOLUTION FINDER ALGORITHMS
% Genetic algorithm
lower = [w_dep(1) w_fb(1) w_arr(1)];           
upper = [w_dep(end) w_fb(end) w_arr(end)];               

% Options for genetic
options_ga = optimoptions('ga', 'PopulationSize', 200, 'FunctionTolerance', 0.01, 'Display', 'off', 'MaxGenerations', 200);

% Solver
N = 19; % Number of departure windows examined
N_ga = 3; % Number of genetic algorithm iteration to have better results
t_opt_ga = [0, 0, 0]; % Storage value for the chosen wind


fprintf('Genetic algorithm computing ... \n\n');
startTime = tic;
[t_opt_ga, dv_min_ga] = ga(@(t) interplanetary1(t(1),t(2),t(3), N_ga), 3, [], [], [], [], lower, upper, [], options_ga);

% Usar el resultado del AG como punto inicial para fmincon
x0 = t_opt_ga; % t_opt es la solución inicial del AG
% Convert dates to Gregorian format
date_dep_ga = mjd20002date(t_opt_ga(1)); % Departure date
date_fb_ga = mjd20002date(t_opt_ga(2)); % Flyby date
date_arr_ga = mjd20002date(t_opt_ga(3)); % Arrival date


%% Refinement with FMINCON // Gradient based optimization 
% fmincon Configuration sqp selection
options_fmincon = optimoptions('fmincon','Display', 'iter-detailed', 'Algorithm', 'sqp','StepTolerance', 1e-10, 'OptimalityTolerance', 1e-6);

fprintf('Refining Solution with FMINCON...\n');
[t_refined, dv_min_refined] = fmincon(@(t) interplanetary1(t(1), t(2), t(3), N_ga),  t_opt_ga, [], [], [], [], lower, upper, [], options_fmincon);

% Convert refined dates to Gregorian format
date_dep_ref = mjd20002date(t_refined(1)); % Departure date
date_fb_ref = mjd20002date(t_refined(2)); % Flyby date
date_arr_ref = mjd20002date(t_refined(3)); % Arrival date


%% Gradient refining method
% Options for gradient
options_grad = optimoptions('fminunc', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals', 1e4, 'MaxIter', 1e4, 'Display', 'off', 'Algorithm', 'quasi-newton'); 

% Solver
[t_opt_grad, dv_min_grad] = fminunc(@(t) interplanetary1(t(1), t(2), t(3), N_ga), t_opt_ga, options_grad);

% Results with gradient
dep_grad = mjd20002date(ceil(t_opt_grad(1)));
fb_grad = mjd20002date(ceil(t_opt_grad(2)));
arr_grad = mjd20002date(ceil(t_opt_grad(3)));

dep_date_grad = sprintf('%02d/%02d/%04d %02d:%02d:%02d', ...
    dep_grad(2), dep_grad(3), dep_grad(1), dep_grad(4), dep_grad(5), dep_grad(6));
fb_date_grad = sprintf('%02d/%02d/%04d %02d:%02d:%02d', ...
    fb_grad(2), fb_grad(3), fb_grad(1), fb_grad(4), fb_grad(5), fb_grad(6));
arr_date_grad = sprintf('%02d/%02d/%04d %02d:%02d:%02d', ...
    arr_grad(2), arr_grad(3), arr_grad(1), arr_grad(4), arr_grad(5), arr_grad(6));

%% Refinamiento Simulated Annealing

options_sa = optimoptions('simulannealbnd', 'MaxIterations', 2000,'Display', 'iter', 'PlotFcns', {@saplotbestx, @saplotbestf, @saplottemperature});

[t_refined_sa, dv_min_sa] = simulannealbnd(@(t) interplanetary1(t(1), t(2), t(3), N_ga), t_opt_ga, lower, upper, options_sa);

date_dep_sa = mjd20002date(t_refined_sa(1)); % Departure date
date_fb_sa = mjd20002date(t_refined_sa(2)); % Flyby date
date_arr_sa = mjd20002date(t_refined_sa(3)); % Arrival date





%% Algoritmo diferencial
options_de = optimoptions('particleswarm','SwarmSize', 200,'MaxIterations', 500,'Display', 'iter', 'UseParallel', true);

% Resolver con PSO / Particle Swarm
[t_opt_pso, dv_min_pso] = particleswarm(@(t) interplanetary1(t(1), t(2), t(3), N_ga),3, lower, upper, options_de);

date_dep_pso = mjd20002date(t_opt_pso(1)); % Departure date
date_fb_pso = mjd20002date(t_opt_pso(2)); % Flyby date
date_arr_pso = mjd20002date(t_opt_pso(3)); % Arrival date

%% Refinement with FMINCON // Gradient based optimization 
% fmincon Configuration sqp selection
options_fmincon = optimoptions('fmincon','Display', 'iter-detailed', 'Algorithm', 'sqp','StepTolerance', 1e-10, 'OptimalityTolerance', 1e-6);

fprintf('Refining Solution with FMINCON...\n');
[t_refined2, dv_min_refined2] = fmincon(@(t) interplanetary1(t(1), t(2), t(3), N_ga),  t_opt_pso, [], [], [], [], lower, upper, [], options_fmincon);

% Convert refined dates to Gregorian format
date_dep_ref2 = mjd20002date(t_refined2(1)); % Departure date
date_fb_ref2 = mjd20002date(t_refined2(2)); % Flyby date
date_arr_ref2 = mjd20002date(t_refined2(3)); % Arrival date


%% Gradient refining method
% Options for gradient
options_grad = optimoptions('fminunc', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals', 1e4, 'MaxIter', 1e4, 'Display', 'off', 'Algorithm', 'quasi-newton'); 

% Solver
[t_opt_grad2, dv_min_grad2] = fminunc(@(t) interplanetary1(t(1), t(2), t(3), N_ga), t_opt_pso, options_grad);

% Results with gradient
dep_grad2 = mjd20002date(ceil(t_opt_grad2(1)));
fb_grad2 = mjd20002date(ceil(t_opt_grad2(2)));
arr_grad2 = mjd20002date(ceil(t_opt_grad2(3)));

dep_date_grad2 = sprintf('%02d/%02d/%04d %02d:%02d:%02d', ...
    dep_grad2(2), dep_grad2(3), dep_grad2(1), dep_grad2(4), dep_grad2(5), dep_grad2(6));
fb_date_grad2 = sprintf('%02d/%02d/%04d %02d:%02d:%02d', ...
    fb_grad2(2), fb_grad2(3), fb_grad2(1), fb_grad2(4), fb_grad2(5), fb_grad2(6));
arr_date_grad2 = sprintf('%02d/%02d/%04d %02d:%02d:%02d', ...
    arr_grad2(2), arr_grad2(3), arr_grad2(1), arr_grad2(4), arr_grad2(5), arr_grad2(6));

%% Refinamiento Simulated Annealing

options_sa = optimoptions('simulannealbnd', 'MaxIterations', 2000,'Display', 'iter', 'PlotFcns', {@saplotbestx, @saplotbestf, @saplottemperature});

[t_refined_sa2, dv_min_sa] = simulannealbnd(@(t) interplanetary1(t(1), t(2), t(3), N_ga), t_opt_pso, lower, upper, options_sa);

date_dep_sa2 = mjd20002date(t_refined_sa2(1)); % Departure date
date_fb_sa2 = mjd20002date(t_refined_sa2(2)); % Flyby date
date_arr_sa2 = mjd20002date(t_refined_sa2(3)); % Arrival date


%% Resultados

% Genetic Algorithm
fprintf('\n Genetic Algorithm Results:\n');
fprintf('Departure: %02d/%02d/%04d\n', date_dep_ga(3), date_dep_ga(2), date_dep_ga(1));
fprintf('Flyby: %02d/%02d/%04d\n', date_fb_ga(3), date_fb_ga(2), date_fb_ga(1));
fprintf('Arrival: %02d/%02d/%04d\n', date_arr_ga(3), date_arr_ga(2), date_arr_ga(1));
fprintf('Delta-v: %.2f km/s\n', dv_min_ga);

% FMINCON Refinement
fprintf('\n FMINCON/Local Refinement Results:\n');
fprintf('Departure: %02d/%02d/%04d\n', date_dep_ref(3), date_dep_ref(2), date_dep_ref(1));
fprintf('Flyby: %02d/%02d/%04d\n', date_fb_ref(3), date_fb_ref(2), date_fb_ref(1));
fprintf('Arrival: %02d/%02d/%04d\n', date_arr_ref(3), date_arr_ref(2), date_arr_ref(1));
fprintf('Delta-v: %.2f km/s\n', dv_min_refined);

% Print refined solution with gradient
fprintf('Refined solution with gradient :\n\n');
fprintf('Departure date : %s \n', dep_date_grad);
fprintf('Fly-by date: %s \n', fb_date_grad);
fprintf('Arrival date : %s \n', arr_date_grad);
fprintf('Minimised cost with gradient : %f km/s \n', dv_min_grad);
fprintf('\n\n')

% Simulated Annealing
fprintf('\n Simulated Annealing Results:\n');
fprintf('Departure: %02d/%02d/%04d \n', date_dep_sa(3), date_dep_sa(2), date_dep_sa(1));
fprintf('Flyby: %02d/%02d/%04d \n', date_fb_sa(3), date_fb_sa(2), date_fb_sa(1));
fprintf('Arrival: %02d/%02d/%04d \n', date_arr_sa(3), date_arr_sa(2), date_arr_sa(1));
fprintf('Delta-v: %.2f km/s\n', dv_min_sa);



% PSO/Differential Algorithm
fprintf('\n PSO/Differential Algorithm Results:\n');
fprintf('Departure: %02d/%02d/%04d\n', date_dep_pso(3), date_dep_pso(2), date_dep_pso(1));
fprintf('Flyby: %02d/%02d/%04d\n', date_fb_pso(3), date_fb_pso(2), date_fb_pso(1));
fprintf('Arrival: %02d/%02d/%04d\n', date_arr_pso(3), date_arr_pso(2), date_arr_pso(1));
fprintf('Delta-v: %.2f km/s\n', dv_min_pso);

% FMINCON Refinement
fprintf('\n FMINCON/Local Refinement Results:\n');
fprintf('Departure: %02d/%02d/%04d\n', date_dep_ref2(3), date_dep_ref2(2), date_dep_ref2(1));
fprintf('Flyby: %02d/%02d/%04d\n', date_fb_ref2(3), date_fb_ref2(2), date_fb_ref2(1));
fprintf('Arrival: %02d/%02d/%04d\n', date_arr_ref2(3), date_arr_ref2(2), date_arr_ref2(1));
fprintf('Delta-v: %.2f km/s\n', dv_min_refined);

% Print refined solution with gradient
fprintf('Refined solution with gradient :\n\n');
fprintf('Departure date : %s \n', dep_date_grad2);
fprintf('Fly-by date: %s \n', fb_date_grad2);
fprintf('Arrival date : %s \n', arr_date_grad2);
fprintf('Minimised cost with gradient : %f km/s \n', dv_min_grad);
fprintf('\n\n')

% Simulated Annealing
fprintf('\n Simulated Annealing Results:\n');
fprintf('Departure: %02d/%02d/%04d \n', date_dep_sa2(3), date_dep_sa2(2), date_dep_sa2(1));
fprintf('Flyby: %02d/%02d/%04d \n', date_fb_sa2(3), date_fb_sa2(2), date_fb_sa2(1));
fprintf('Arrival: %02d/%02d/%04d \n', date_arr_sa2(3), date_arr_sa2(2), date_arr_sa2(1));
fprintf('Delta-v: %.2f km/s\n', dv_min_sa);

%% GENERATE CHOP PLOT
fprintf('Generating CHOP Plot...\n');

% Convert MJD2000 to Gregorian dates for departure and arrival
t_dep_greg = arrayfun(@mjd20002date, w_dep, 'UniformOutput', false); % Convert each MJD2000 in t_dep
t_arr_greg = arrayfun(@mjd20002date, w_arr, 'UniformOutput', false); % Convert each MJD2000 in t_arr

% Convert the cell arrays to readable format for plotting
t_dep_labels = cellfun(@(x) datestr(datenum(x(1:3))), t_dep_greg, 'UniformOutput', false); % Convert to date strings
t_arr_labels = cellfun(@(x) datestr(datenum(x(1:3))), t_arr_greg, 'UniformOutput', false); % Convert to date strings


% Initialize delta-v matrix
delta_v = NaN(length(w_dep), length(w_arr));

% Calculate delta-v for all combinations of departure and arrival dates
for i = 1:length(w_dep)
    for j = 1:length(w_arr)
        if w_arr(j) > w_dep(i) % Ensure arrival is after departure
            % Call the interplanetary function
            delta_v(i, j) = interplanetary1(w_dep(i), (w_dep(i) + w_arr(j)) / 2, w_arr(j),N_ga);
        end
    end
end


% Plot the CHOP Plot
figure;
hold on;
contourf(w_dep, w_arr, delta_v', 20, 'LineColor', 'none'); % Contour plot
colorbar; % Add a color bar
colormap('jet'); % Use jet colormap

% Customize the plot
xlabel('Departure Time');
ylabel('Arrival Time');
title('Delta-v Contour Plot');

% Update X and Y ticks with readable labels
set(gca, 'XTick', w_dep(1:30:end), 'XTickLabel', t_dep_labels(1:30:end));
set(gca, 'YTick', w_arr(1:30:end), 'YTickLabel', t_arr_labels(1:30:end));

% Highlight a specific chosen window
chosen_departure = w_dep(round(length(w_dep) / 2)); % Example chosen departure date
chosen_arrival = w_arr(round(length(w_arr) / 2));  % Example chosen arrival date
plot(chosen_departure, chosen_arrival, 'p', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'cyan'); % Add marker
legend('Chosen time window', 'Location', 'northeast');

hold off;

%% PLOT RESULTS
%% Results 
[dv_opt, dv_dep, dv_arr, r1, v1i, r2, v2f, r3, v3f, v1t, v2t, v2t_1, v3t, vinfmin_vec, vinfplus_vec] = interplanetary1(t_opt_grad(1), t_opt_grad(2), t_opt_grad(3),N_ga);
[vinfm, vinfp, delta, rp, am, ap, em, ep, vpm, vpp, deltam, deltap, dv_fb_tot, dv_fb_pow] = flyby_powered(vinfmin_vec, vinfplus_vec, muE);

%fprintf('The final solution is :\n\n');
%fprintf('Departure date : %s \n', dep_date_grad);
%fprintf('Fly-by date: %s \n', fb_date_grad);
%fprintf('Arrival date : %s \n', arr_date_grad);
%fprintf('\n\n');

%% Heliocentric trajectory
% Initialisation
N_t = 50000;

t_dep = t_opt_grad(1) * 86400;
t_fb = t_opt_grad(2) * 86400;
t_arr = t_opt_grad(3) * 86400;

dt_leg1 = t_fb - t_dep;
dt_leg2 = t_arr - t_fb;

tspan_mercury = linspace(0, T_dep, N_t);
tspan_leg1 = linspace(0, -dt_leg1, N_t);
tspan_earth = linspace(0, T_fb, N_t);
tspan_leg2 = linspace(0, -dt_leg2, N_t);
tspan_asteroid = linspace(0, T_arr, N_t);

% Set options for ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Matrices defining
y_mercury = [ r1; v1i ];
y_leg1 = [ r2; v2t' ];
y_earth = [ r2; v2f ];
y_leg2 = [ r3; v3t' ];
y_ast = [ r3; v3f ];

[ t1, Y_mercury ] = ode113(@(t,y) ode_2bp(t,y,muS), tspan_mercury, y_mercury, options);
[ t2, Y_leg1 ] = ode113(@(t,y) ode_2bp(t,y,muS), tspan_leg1, y_leg1, options);
[ t3, Y_earth ] = ode113(@(t,y) ode_2bp(t,y,muS), tspan_earth, y_earth, options);
[ t4, Y_leg2 ] = ode113(@(t,y) ode_2bp(t,y,muS), tspan_leg2, y_leg2, options);
[ t5, Y_ast] = ode113(@(t,y) ode_2bp(t,y,muS), tspan_asteroid, y_ast, options);

% Plot
n = 1.496e+8;

figure();
plot3(Y_mercury(:,1)/n, Y_mercury(:,2)/n,  Y_mercury(:,3)/n, 'b-', 'LineWidth', 1);
hold on;
plot3(Y_leg1(:,1)/n, Y_leg1(:,2)/n,  Y_leg1(:,3)/n, 'm-', 'LineWidth', 1);
plot3(Y_earth(:,1)/n, Y_earth(:,2)/n,  Y_earth(:,3)/n, 'r-', 'LineWidth', 1);
plot3(Y_leg2(:,1)/n, Y_leg2(:,2)/n,  Y_leg2(:,3)/n, 'g-', 'LineWidth', 1);
plot3(Y_ast(:,1)/n, Y_ast(:,2)/n, Y_ast(:,3)/n, 'y-', 'LineWidth', 1);
xlabel('X [AU]');
ylabel('Y [AU]');
zlabel('Z [AU]');
title('Heliocentric trajectory');
axis([-2.5 2.5 -2.5 2.5 -2.5 2.5]);
grid on;

legend("Mercury's orbit", "Transfer orbit to Earth", "Earth's orbit", "Transfer orbit to the asteroid", "Asteroid's orbit");
hold off;

%% Fly-by trajectory (planetocentric)
% Results
IPm = -am*em*cos(deltam); % Impact paramater for incoming hyperbolic trajectory
IPp = -ap*ep*cos(deltap); % Impact paramater for outgoing hyperbolic trajectory

CA = min(IPm, IPp); % Altitude of the closest approach

fprintf('The altitude of the closest approach is : %f km \n\n', CA);

fprintf('The total velocity change due to flyby is : %f km/s \n', dv_fb_tot);
fprintf('The cost of the manoeuvre at pericentre is : %f km/s \n\n', dv_fb_pow);

% Initial conditions planetocentric
u = cross(vinfmin_vec,vinfplus_vec)/norm(cross(vinfmin_vec,vinfplus_vec));

betam = pi/2 - deltam/2;

dir_vm = vinfmin_vec/norm(vinfmin_vec); % Vinf- velocity direction
dir_vp = vinfplus_vec/norm(vinfplus_vec); % Vinf+ velocity direction

dirm = Rotate(dir_vm, u, deltam/2); 
dirp = Rotate(dir_vp, u, -deltap/2);

r0 = rp * Rotate(dir_vm, u, -betam);

vm = vpm*dirm;
vp = vpp*dirp;

% Time span planetocentric
tspan_m = linspace(0, -50000, 100000);
tspan_p = linspace(0, 50000, 100000);

% Set options for ODE solver
options_fb = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

% Integration of planetocentric trajectory
y0m = [r0; vm];
y0p = [r0; vp];

[t_fb_min, Y_fb_min] = ode113(@(t, y) ode_2bp(t, y, muE), tspan_m, y0m, options_fb);

[t_fb_plus, Y_fb_plus] = ode113(@(t, y) ode_2bp(t, y, muE), tspan_p, y0p, options_fb);

% Plot
figure();
hold on;

plot(Y_fb_min(:, 1) / Re, Y_fb_min(:, 2) / Re, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Flyby hyperbola (infront)');
plot(Y_fb_plus(:, 1) / Re, Y_fb_plus(:, 2) / Re, 'b-', 'LineWidth', 1.5);
plot(0, 0, 'yo', 'MarkerSize', 15, 'MarkerFaceColor', 'blue', 'DisplayName', 'Earth');

xlabel('x [Re]');
ylabel('y [Re]');
title('Trajectory in Earth-centred frame parallel to (HECI)');
axis equal;
grid on;

xlim([-6, 9]);
ylim([-9, 3]);