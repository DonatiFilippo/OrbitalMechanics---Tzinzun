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
% LAST UPDATE: 27-12-2024
%
%--------------------------------------------------------------------------

%% SCRIPT INITIALISATION
clc
clear variables
close all

%% Saving data
if exist(fullfile(cd, 'Assigment1.cvs'), 'file') == 2
    delete(fullfile(cd, 'Assigment1.txt'))
end

filename = 'Assigment1.txt';
fileID = fopen(filename,'w+');
fprintf(fileID,'[TEAM 2427 - TZINTZUN]\n');
fclose(fileID);

%% PRESENTATION
fileID = fopen(filename,'a+');
fprintf(fileID,'----------------------------------------------\n');
fprintf(fileID,'     INTERPLANETARY MISSION ASSIGNMENT        \n');
fprintf(fileID,'----------------------------------------------\n');
fprintf(fileID,'\n');
fclose(fileID);

%% IMPOSED DATES
Departure_date = [2030, 1, 1, 0, 0, 0]; % Earliest departure date [Gregorian date]
Arrival_date = [2060, 1, 1, 0, 0, 0];   % Latest arrival date [Gregorian date]
 
mjd_dep = date2mjd2000(Departure_date); % Earliest departure date converted to modified Julian day 2000
mjd_arr = date2mjd2000(Arrival_date);   % Latest arrival date converted to modified Julian day 2000

fileID = fopen(filename,'a+');
fprintf(fileID,'Mercury to Asteroid no.30 w flyby on Earth\n');
fprintf(fileID,'Mission Window : [%d %d %d %d %d %d] - [%d %d %d %d %d %d]\n', Departure_date, Arrival_date);
fclose(fileID);

%% PHYSICAL PARAMETERS
Departure_planet = 1;     % Mercury as the departure planet 
Flyby_planet = 3;         % Earth as the flyby planet
Arrival_asteroid_id = 30; % Asteroid no.30 as the arrival objective

muM = astroConstants(11);  % Mercury's gravitational constant [km^3/s^2]
muE = astroConstants(13);  % Earth's gravitational constant [km^3/s^2]
muS = astroConstants(4);   % Sun's gravitational constant [km^3/s^2]
Re  = astroConstants(23);  % Earth's mean radius [km]

% Keplerian elements computing / kep = [a e i Om om theta] [km, rad]
[kep_dep, ~] = uplanet(0, Departure_planet);            % Mercury's keplerian elements at initial time
[kep_fb, ~] = uplanet(0, Flyby_planet);                 % Earth's keplerian elements at initial time
[kep_arr, ~, ~] = ephAsteroids(0, Arrival_asteroid_id); % Asteroid's keplerian elements at initial time

%% PERIODS COMPUTING
T_dep = 2*pi*sqrt(kep_dep(1)^3/muS); % Mercury's orbital period [s]
T_fb  = 2*pi*sqrt(kep_fb(1)^3/muS);  % Earth's orbital period [s]
T_arr = 2*pi*sqrt(kep_arr(1)^3/muS); % Asteroid's orbital period [s]

T_syn_dep2fb  = T_dep * T_fb/abs(T_dep - T_fb);   % Mercury's synodic orbital period with respect to Earth [s]
T_syn_fb2arr  = T_fb * T_arr/abs(T_fb - T_arr);   % Earth's synodic orbital period with respect to the asteroid [s]
T_syn_dep2arr = T_dep * T_arr/abs(T_dep - T_arr); % Mercury's synodic orbital period with respect to the asteroid [s]

fileID = fopen(filename,'a+');
fprintf(fileID,'Displaying synodics period in years :\n\n');
fprintf(fileID,'Mercury''s synodic orbital period with respect to Earth : %f years \n', T_syn_dep2fb/(86400*365.25));
fprintf(fileID,'Earth''s synodic orbital period with respect to the asteroid : %f years \n', T_syn_fb2arr/(86400*365.25));
fprintf(fileID,'Mercury''s synodic orbital period with respect to the asteroid : %f years \n', T_syn_dep2arr/(86400*365.25));
fprintf(fileID,'\n\n');
fclose(fileID);

%% WINDOWS RESEARCH

a_t1 = (kep_dep(1) + kep_fb(1))/2; % Semi-major axis of the first transfert arc between Mercury and Earth [km]
a_t2 = (kep_fb(1) + kep_arr(1))/2; % Semi-major axis of the second transfert arc between Earth and the asteroid [km]

T_t1 = 2*pi*sqrt(a_t1^3/muS); % Period of the first transfer arc [s]
T_t2 = 2*pi*sqrt(a_t2^3/muS); % Period of the second transfer arc [s]

tof_t1 = 1/2*T_t1 / 86400;  % Time of flight of the first transfer arc [days]
tof_t2 = 1/2*T_t2 / 86400;  % Time of flight of the second transfer arc [days]
tof_h  = (tof_t1 + tof_t2); % Total time of flight for Hohmann transfer [days]

fileID = fopen(filename,'a+');
fprintf(fileID,'Displaying tof of transfer arc for Hohmann transfers in years :\n\n');
fprintf(fileID,'Time of flight of the first transfer arc from Mercury to Earth : %f years\n', tof_t1/(365.25));
fprintf(fileID,'Time of flight of the second transfer arc from Earth to the asteroid : %f years \n', tof_t2/(365.25));
fprintf(fileID,'Total time of flight from Mercury to the asteroid : %f years \n', tof_h/(365.25));
fprintf(fileID,'\n\n');
fclose(fileID);

step = 1;   % Step size of days for iterating through time windows
SM   = 0.4; % Safety margin for time of flight, 40% adjustment based on bibliography

% Calculate the synodic period with the most relevance
SP = max([T_syn_dep2fb, T_syn_fb2arr, T_syn_dep2arr]) / 86400; % Synodic period in days

% Time of flight ranges with safety margins
tof_t1_min = (1 - SM) * tof_t1; % Minimum time of flight Mercury -> Earth
tof_t1_max = (1 + SM) * tof_t1; % Maximum time of flight Mercury -> Earth

tof_t2_min = (1 - SM) * tof_t2; % Minimum time of flight Earth -> Asteroid
tof_t2_max = (1 + SM) * tof_t2; % Maximum time of flight Earth -> Asteroid

% Calculate the last possible departure time from Mercury
t_ldM = SP - tof_t1_min; % Last departure from Mercury to arrive within the synodic period

% Define the departure window from Mercury
w_dep = mjd_dep : step : mjd_dep + t_ldM; % First departure window from Mercury

% Calculate the arrival window at Earth
w_fb_min = w_dep(1) + tof_t1_min;   % Earliest arrival at Earth
w_fb_max = w_dep(end) + tof_t1_max; % Latest arrival at Earth
w_fb = w_fb_min : step : w_fb_max;  % Arrival window at Earth

% Calculate the departure window from Earth to the asteroid
w_arr_min = w_fb(1) + tof_t2_min;     % Earliest departure from Earth to the asteroid
w_arr_max = w_fb(end) + tof_t2_max;   % Ensure compatibility with the synodic period
w_arr = w_arr_min : step : w_arr_max; % Arrival window at the asteroid

%% BEST SOLUTION FINDER ALGORITHMS
%% Genetic algorithm
lower = [w_dep(1) w_fb(1) w_arr(1)];           
upper = [w_dep(end) w_fb(end) w_arr(end)];               

lower_ga = [w_dep(1) w_fb(1) w_arr(1)];
upper_ga = [w_dep(end) w_fb(end) w_arr(end)];

% Options for genetic
options_ga = optimoptions('ga', 'PopulationSize', 300, ...
    'FunctionTolerance', 0.001, 'Display', 'off', 'MaxGenerations', 200);
 
% Solver for genetic
N         = ceil((mjd_arr-w_arr_max)/365.25);
N_ga      = 3;         % Number of genetic algorithm iteration to have better results
dv_min_ga = 50;        % Arbitrary chosen value of total cost
t_opt_ga  = [0, 0, 0]; % Storage value for the chosen windows

fprintf('Genetic algorithm computing ... \n\n');
startTime = tic;
for i = 1:N
    fprintf('ITERATION NUMBER : %2.f \n \n', i);

    for j = 1:N_ga
        [t_opt_ga_computed, dv_min_ga_computed] = ga(@(t) interplanetary(t(1),t(2),t(3)), 3, [], [], [], [], lower, upper, [], options_ga);
        if dv_min_ga_computed < dv_min_ga && t_opt_ga_computed(3) < mjd_arr
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
date_dep_ga = mjd20002date(t_opt_ga(1));
date_fb_ga  = mjd20002date(t_opt_ga(2));
date_arr_ga = mjd20002date(t_opt_ga(3));

%% Refinement with FMINCON
% fmincon Configuration sqp selection options
options_fmincon = optimoptions('fmincon','Display', 'off', 'Algorithm', 'sqp','StepTolerance', 1e-10, 'OptimalityTolerance', 1e-8);

% Fmincon solver
fprintf('Refining Solution with FMINCON...\n');
[t_refined_fmin, dv_min_fmin] = fmincon(@(t) interplanetary(t(1), t(2), t(3)),  t_opt_ga, [], [], [], [], lower_ga, upper_ga, [], options_fmincon);

% Convert refined dates to Gregorian format
date_dep_ref = mjd20002date(t_refined_fmin(1));
date_fb_ref  = mjd20002date(t_refined_fmin(2));
date_arr_ref = mjd20002date(t_refined_fmin(3));

%% Gradient refining method
% Options for gradient
options_grad = optimoptions('fminunc', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals', 1e4, 'MaxIter', 1e4, 'Display', 'off', 'Algorithm', 'quasi-newton'); 

% Gradient solver
[t_refined_grad, dv_min_grad] = fminunc(@(t) interplanetary(t(1), t(2), t(3)), t_opt_ga, options_grad);

% Results with gradient
date_dep_grad = mjd20002date(t_refined_grad(1));
date_fb_grad  = mjd20002date(t_refined_grad(2));
date_arr_grad = mjd20002date(t_refined_grad(3));

%% Simulated annealing refining method
% Options for simulated annealing
options_sa = optimoptions('simulannealbnd', 'MaxIterations', 2000, 'Display', 'off', 'PlotFcns', []);

% Simulated annealing solver
[t_refined_sa, dv_min_sa] = simulannealbnd(@(t) interplanetary(t(1), t(2), t(3)), t_opt_ga, lower_ga, upper_ga, options_sa);

% Results with simulated annealing
date_dep_sa = mjd20002date(t_refined_sa(1));
date_fb_sa = mjd20002date(t_refined_sa(2));
date_arr_sa = mjd20002date(t_refined_sa(3));

%% Algorithm comparison
% Genetic Algorithm
fileID = fopen(filename,'a+');
fprintf(fileID,'\n\n')
fprintf(fileID,'Genetic Algorithm Results:\n\n');
fprintf(fileID,'Departure: %02d/%02d/%04d\n', date_dep_ga(3), date_dep_ga(2), date_dep_ga(1));
fprintf(fileID,'Flyby: %02d/%02d/%04d\n', date_fb_ga(3), date_fb_ga(2), date_fb_ga(1));
fprintf(fileID,'Arrival: %02d/%02d/%04d\n', date_arr_ga(3), date_arr_ga(2), date_arr_ga(1));
fprintf(fileID,'Delta-v: %.2f km/s\n', dv_min_ga);
fprintf(fileID,'\n\n');
fclose(fileID);

% FMINCON Refinement
fileID = fopen(filename,'a+');
fprintf(fileID,'FMINCON/Local Refinement Results:\n\n');
fprintf(fileID,'Departure: %02d/%02d/%04d\n', date_dep_ref(3), date_dep_ref(2), date_dep_ref(1));
fprintf(fileID,'Flyby: %02d/%02d/%04d\n', date_fb_ref(3), date_fb_ref(2), date_fb_ref(1));
fprintf(fileID,'Arrival: %02d/%02d/%04d\n', date_arr_ref(3), date_arr_ref(2), date_arr_ref(1));
fprintf(fileID,'Delta-v: %.2f km/s\n', dv_min_fmin);
fprintf(fileID,'\n\n');
fclose(fileID);

% Refined solution with gradient
fileID = fopen(filename,'a+');
fprintf(fileID,'Refined solution with gradient :\n\n');
fprintf(fileID,'Departure date : %02d/%02d/%04d\n', date_dep_grad(3), date_dep_grad(2), date_dep_grad(1));
fprintf(fileID,'Fly-by date: %02d/%02d/%04d\n', date_fb_grad(3), date_fb_grad(2), date_fb_grad(1));
fprintf(fileID,'Arrival date : %02d/%02d/%04d\n', date_arr_grad(3), date_arr_grad(2), date_arr_grad(1));
fprintf(fileID,'Minimised cost with gradient : %f km/s \n', dv_min_grad);
fprintf(fileID,'\n\n');
fclose(fileID);

% Simulated Annealing
fileID = fopen(filename,'a+');
fprintf(fileID,'Simulated Annealing Results:\n\n');
fprintf(fileID,'Departure: %02d/%02d/%04d \n', date_dep_sa(3), date_dep_sa(2), date_dep_sa(1));
fprintf(fileID,'Flyby: %02d/%02d/%04d \n', date_fb_sa(3), date_fb_sa(2), date_fb_sa(1));
fprintf(fileID,'Arrival: %02d/%02d/%04d \n', date_arr_sa(3), date_arr_sa(2), date_arr_sa(1));
fprintf(fileID,'Delta-v: %.2f km/s\n', dv_min_sa);
fprintf(fileID,'\n\n');
fclose(fileID);

%% Choice of the best solution

dv_min_sol = min([dv_min_ga, dv_min_fmin, dv_min_sa, dv_min_grad]);

if dv_min_sol == dv_min_fmin
    t_opt_sol = t_refined_fmin;
elseif dv_min_sol == dv_min_grad
    t_opt_sol = t_refined_grad;
elseif dv_min_sol == dv_min_sa
    t_opt_sol = t_refined_sa;
end

date_dep_sol = mjd20002date(t_opt_sol(1));
date_fb_sol  = mjd20002date(t_opt_sol(2));
date_arr_sol = mjd20002date(t_opt_sol(3));

%% PLOT RESULTS
%% Results 
[dv_opt, dv_dep, dv_arr, r1, v1i, r2, v2f, r3, v3f, v1t, v2t, v2t_1, v3t, vinfmin_vec, vinfplus_vec, ~] = interplanetary(t_opt_sol(1), t_opt_sol(2), t_opt_sol(3));
[vinfm, vinfp, delta, rp, am, ap, em, ep, vpm, vpp, deltam, deltap, dv_fb_tot, dv_fb_pow]               = flyby_powered(vinfmin_vec, vinfplus_vec, muE);

v1_t = v1t';        
v2_t = v2t_1';  

[a1, e1, i1, Omega1, omega1, nu1] = car2kep(r1, v1_t);
[a2, e2, i2, Omega2, omega2, nu2] = car2kep(r2, v2_t); 

fileID = fopen(filename,'a+');
fprintf(fileID,'The final solutions are :\n\n');
fprintf(fileID,'Departure: %02d/%02d/%04d - %02d:%02d:%02d \n', date_dep_sol(3), date_dep_sol(2), date_dep_sol(1), date_dep_sol(4), date_dep_sol(5),date_dep_sol(6));
fprintf(fileID,'Flyby: %02d/%02d/%04d - %02d:%02d:%02d \n', date_fb_sol(3), date_fb_sol(2), date_fb_sol(1), date_fb_sol(4), date_fb_sol(5), date_fb_sol(6));
fprintf(fileID,'Arrival: %02d/%02d/%04d - %02d:%02d:%02d \n', date_arr_sol(3), date_arr_sol(2), date_arr_sol(1), date_arr_sol(4), date_arr_sol(5), date_arr_sol(6));
fprintf(fileID,'Mercury to Earth Orbit: ');
fprintf(fileID,'a = %.2f km, e = %.5f, i = %.2f°, Ω = %.2f°, ω = %.2f°, ν = %.2f°\n', a1, e1, i1, Omega1, omega1, nu1);
fprintf(fileID,'Earth to Asteroid Orbit: ');
fprintf(fileID,'a = %.2f km, e = %.5f, i = %.2f°, Ω = %.2f°, ω = %.2f°, ν = %.2f°\n', a2, e2, i2, Omega2, omega2, nu2);
fprintf(fileID,'Delta-v at departure: %.2f km/s\n', dv_dep);
fprintf(fileID,'Delta-v at arrival: %.2f km/s\n', dv_arr);
fprintf(fileID,'Delta-v optimal: %.2f km/s\n', dv_opt);
fprintf(fileID,'Delta-v minimun: %.2f km/s\n', dv_min_sol);
fprintf(fileID,'\n\n');
fclose(fileID);

%% Heliocentric trajectory
% Initialisation
N_t = 50000;

t_dep = t_refined_grad(1) * 86400;
t_fb  = t_refined_grad(2) * 86400;
t_arr = t_refined_grad(3) * 86400;

dt_leg1 = t_fb - t_dep;
dt_leg2 = t_arr - t_fb;

tspan_mercury  = linspace(0, T_dep, N_t);
tspan_leg1     = linspace(0, -dt_leg1, N_t);
tspan_earth    = linspace(0, T_fb, N_t);
tspan_leg2     = linspace(0, -dt_leg2, N_t);
tspan_asteroid = linspace(0, T_arr, N_t);

% Set options for ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Matrices defining
y_mercury = [ r1; v1i ];
y_leg1    = [ r2; v2t'];
y_earth   = [ r2; v2f ];
y_leg2    = [ r3; v3t'];
y_ast     = [ r3; v3f ];

[ t1, Y_mercury ] = ode113(@(t,y) ode_2bp(t,y,muS), tspan_mercury, y_mercury, options);
[ t2, Y_leg1    ] = ode113(@(t,y) ode_2bp(t,y,muS), tspan_leg1, y_leg1, options);
[ t3, Y_earth   ] = ode113(@(t,y) ode_2bp(t,y,muS), tspan_earth, y_earth, options);
[ t4, Y_leg2    ] = ode113(@(t,y) ode_2bp(t,y,muS), tspan_leg2, y_leg2, options);
[ t5, Y_ast     ] = ode113(@(t,y) ode_2bp(t,y,muS), tspan_asteroid, y_ast, options);

% Plot
n = astroConstants(2);

figure('Name', 'Heliocentric trajectory', 'NumberTitle', 'on', 'Position', [0, 250, 500, 500], 'Color', [1 1 1]);
plot3(Y_mercury(:,1)/n, Y_mercury(:,2)/n, Y_mercury(:,3)/n, 'LineWidth', 1, 'Color', [0.4039 0.5333 0.6353]);
hold on;
plot3(Y_leg1(:,1)/n, Y_leg1(:,2)/n, Y_leg1(:,3)/n,'LineWidth', 1, 'Color', [0.9290 0.6940 0.1250]); 
plot3(Y_earth(:,1)/n, Y_earth(:,2)/n, Y_earth(:,3)/n,'LineWidth', 1, 'Color', [0.4980 0.6706 0.5255]); 
plot3(Y_leg2(:,1)/n, Y_leg2(:,2)/n, Y_leg2(:,3)/n, 'LineWidth', 1, 'Color', [0.4940 0.1840 0.5560]); 
plot3(Y_ast(:,1)/n, Y_ast(:,2)/n, Y_ast(:,3)/n, 'LineWidth', 1, 'Color', [0.8667 0.5608 0.4314]); 


Planet3d(10, [0, 0, 0], '~', 0.0931);                                                       % Sun scaled to 20 times its size in AU for visualization purposes only.
Planet3d(2, [Y_mercury(end,1)/n, Y_mercury(end,2)/n, Y_mercury(24833,3)/n], '~', 0.0489);   % Mercury scaled to 3000 times its size in AU for visualization purposes only.
Planet3d(0, [Y_earth(500,1)/n, Y_earth(500,2)/n, Y_earth(1,3)/n], '~', 0.0426);             % Eart scaled to 1000 times its size in AU for visualization purposes only.
Planet3d(11, [Y_ast(end,1)/n, Y_ast(end,2)/n, Y_ast(1,3)/n], '~', 0.0200);                  % Asteroid just for visualization purposes only.


xlabel('X [AU]');
ylabel('Y [AU]');
zlabel('Z [AU]');
title('Heliocentric trajectory');
axis([-2.5 2.5 -2.5 2.5 -2.5 2.5]);
grid on;

legend("Mercury's orbit", "Transfer orbit to Earth", "Earth's orbit", "Transfer orbit to the asteroid", "Asteroid's orbit", 'Location', 'northeast');
hold off;

%% Fly-by trajectory (planetocentric)
% Results
CA = rp - Re; % Altitude of closest approach

fileID = fopen(filename,'a+');
fprintf(fileID,'The altitude of the closest approach is : %f km \n\n', CA);
fprintf(fileID,'The total velocity change due to flyby is : %f km/s \n', dv_fb_tot);
fprintf(fileID,'The cost of the manoeuvre at pericentre is : %f km/s \n\n', dv_fb_pow);
fclose(fileID);

% Initial conditions planetocentric
u = cross(vinfmin_vec,vinfplus_vec)/norm(cross(vinfmin_vec,vinfplus_vec));

betam = pi/2 - deltam/2;

dir_vm = vinfmin_vec/norm(vinfmin_vec);   % Vinf- velocity direction
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

[t_fb_min, Y_fb_min  ] = ode113(@(t, y) ode_2bp(t, y, muE), tspan_m, y0m, options_fb);
[t_fb_plus, Y_fb_plus] = ode113(@(t, y) ode_2bp(t, y, muE), tspan_p, y0p, options_fb);

% Plot
figure('Name', 'Fly-by trajectory (planetocentric)', 'NumberTitle', 'on', 'Position', [500, 250, 400, 400], 'Color', [1 1 1]);
hold on;

plot3(Y_fb_min(:, 1) / Re, Y_fb_min(:, 2) / Re, Y_fb_min(:, 3) / Re, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5, 'DisplayName', 'Flyby hyperbola (infront)');
plot3(Y_fb_plus(:, 1) / Re, Y_fb_plus(:, 2) / Re, Y_fb_plus(:, 3) / Re, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.5);
Planet3d(0, [0, 0, 0], 'RE'); 

view(3);

xlabel('x [Re]');
ylabel('y [Re]');
zlabel('z [Re]');
title('Trajectory in Earth-centred frame parallel to (HECI)');
axis equal;
grid on;

xlim([-10, 10]);
ylim([-10, 10]);
zlim([-10, 10]);

%% Porkchop patches
% Windows for pork-chop to work
step_plot = 25;

w_dep_plot = lower_ga(1) : step_plot : upper_ga(1);
 w_fb_plot = lower_ga(2) : step_plot : upper_ga(2);
w_arr_plot = lower_ga(3) : step_plot : upper_ga(3);

[delta_v1, delta_v2] = porkchop(w_dep_plot, w_fb_plot, w_fb_plot, w_arr_plot);

[w_dep_datenum, w_fb_datenum, w_arr_datenum, t_opt_datenum] = datedata(w_dep_plot, w_fb_plot, w_arr_plot, t_opt_sol);

figure('Name', 'Pork chop plot contour from Mercury to Earth', 'NumberTitle', 'on', 'Position', [900, 0, 400, 350], 'Color', [1 1 1])
hold on
grid on
title('Pork chop plot contour from Mercury to Earth')
xlabel('Time of arrival [MJD2000]')
ylabel('Time of departure [MJD2000]')
contour(w_fb_datenum, w_dep_datenum, delta_v1, 50);
colorbar
colormap jet
datetick('x', 'dd/mm/yyyy', 'keepticks', 'keeplimits')
datetick('y', 'dd/mm/yyyy', 'keepticks', 'keeplimits')
set(gca, 'XTickLabelRotation', 45)
set(gca, 'YTickLabelRotation', 45)
%hold on
plot(t_opt_datenum(1), t_opt_datenum(2), 'ro', 'MarkerSize', 10)
legend('Chosen time window');
hold off

figure('Name', 'Pork chop plot contour from Earth to Asteroid', 'NumberTitle', 'on', 'Position', [900, 300, 400, 350], 'Color', [1 1 1])
hold on
grid on
title('Pork chop plot contour from Earth to Asteroid')
xlabel('Time of arrival [MJD2000]')
ylabel('Time of departure [MJD2000]')
contour(w_arr_datenum, w_fb_datenum, delta_v2, 50);
colorbar
colormap jet
datetick('x', 'dd/mm/yyyy', 'keepticks', 'keeplimits')
datetick('y', 'dd/mm/yyyy', 'keepticks', 'keeplimits')
set(gca, 'XTickLabelRotation', 45)
set(gca, 'YTickLabelRotation', 45)
plot(t_opt_datenum(3), t_opt_datenum(2), 'ro', 'MarkerSize', 10)
legend('Chosen time window');
hold off

%% PorkChop plot

% 3D Pork Chop Plot: Mercury to Earth
figure('Name', '3D visualization of Transfer Delta-V: Mercury to Earth', 'NumberTitle', 'on', 'Position', [0, -30, 500, 300], 'Color', [1 1 1])
hold on
grid on
title('3D Pork chop plot: Mercury to Earth')
xlabel('Time of arrival [MJD2000]')
ylabel('Time of departure [MJD2000]')
zlabel('Delta V [km/s]')
[X, Y] = meshgrid(linspace(min(w_arr_datenum), max(w_arr_datenum), 50), ...
                  linspace(min(w_fb_datenum), max(w_fb_datenum), 50));
delta_v_interp = interp2(w_arr_datenum, w_fb_datenum, delta_v1, X, Y, 'spline');
surf(X, Y, delta_v_interp, 'EdgeColor', 'k', 'FaceAlpha', 0.5);
colormap parula
colorbar
%caxis([min(delta_v1(:)), max(delta_v1(:))])
datetick('x', 'dd/mm/yyyy', 'keepticks', 'keeplimits')
datetick('y', 'dd/mm/yyyy', 'keepticks', 'keeplimits')
set(gca, 'XTickLabelRotation', 45)
set(gca, 'YTickLabelRotation', 45)
%plot3(t_opt_datenum(1), t_opt_datenum(2), min(delta_v1(:)), 'ro', 'MarkerSize', 10, 'LineWidth', 2)
view(3)
%legend('Chosen time window', 'Location', 'northeast');
hold off

 
figure('Name', '3D visualization of Transfer Delta-V: Earth to Asteroid', 'NumberTitle', 'on', 'Position', [400, -30, 500, 300], 'Color', [1 1 1])
hold on
grid on
title('3D Pork chop plot: Earth to Asteroid')
xlabel('Time of arrival [MJD2000]')
ylabel('Time of departure [MJD2000]')
zlabel('Delta V [km/s]')
[X, Y] = meshgrid(linspace(min(w_arr_datenum), max(w_arr_datenum), 50), ...
                  linspace(min(w_fb_datenum), max(w_fb_datenum), 50));
delta_v_interp = interp2(w_arr_datenum, w_fb_datenum, delta_v2, X, Y, 'spline');
surf(X, Y, delta_v_interp, 'EdgeColor', 'k', 'FaceAlpha', 0.5);
colormap parula
colorbar
%caxis([min(delta_v2(:)), max(delta_v2(:))])
datetick('x', 'dd/mm/yyyy', 'keepticks', 'keeplimits')
datetick('y', 'dd/mm/yyyy', 'keepticks', 'keeplimits')
set(gca, 'XTickLabelRotation', 45)
set(gca, 'YTickLabelRotation', 45)
%plot3(t_opt_datenum(3), t_opt_datenum(2), min(delta_v2(:)), 'ro', 'MarkerSize', 10, 'LineWidth', 2)
view(3)
%legend('Chosen time window', 'Location', 'northeast');
hold off
