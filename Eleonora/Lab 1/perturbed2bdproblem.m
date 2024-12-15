% PERTURBATED TWO-BODY PROBLEM
clear;
clc;
close all;

% Parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
Rt = astroConstants(23); % Earth's radius
J2 = astroConstants(9); % second zonal harmonic

% Initial conditions
% HIGHLY ECCENTRIC AND INCLINED ORBIT
r0 = [6495; -970; -3622]; % [km]
v0 = [4.752; 2.130; 7.950]; % [km/s]
y0 = [r0; v0]; % [6x1]

% QUASI-CIRCULAR MEDIUM ORBIT
% r0 = [26578.137; 0; 0]; % [km]
% v0 = [0; 2.221; 3.173];  % [km/s]

% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [s]
tspan = linspace( 0, 730*T, 100000 );
l = length(tspan);

% Options for the ODE solver
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

% Perform the integration for the perturbed case
[Tout_p, Yout_p] = ode113 (@(t,y) ode_2bp_pert(t,y,mu_E,Rt,J2), tspan, y0, options);
Tout_p = Tout_p'; % [1xl] row vector
Yout_p = Yout_p'; % [6xl] row vector

% Perform the integration for the standard (non perturbed) case
[Tout, Yout] = ode113 (@(t,y) ode_2bp(t,y,mu_E), tspan, y0, options);
Tout = Tout'; % [1xl] row vector
Yout = Yout'; % [6xl] row vector

% Plot the results
Earth_3D(Rt)
hold on
plot3( Yout(1, :), Yout(2, :), Yout(3, :), 'b-', 'LineWidth', 2)
plot3( Yout_p(1, :), Yout_p(2, :), Yout_p(3, :),'g-')
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;

%% Angular momentum and eccentricity plots

% We want to check that h and e remain constant in magnitude and direction 

% Dimensions definitions
h_vect_p = zeros(3, l); 
h_p = zeros(1, l); 
e_vect_p = zeros(3, l);
e_p = zeros(1, l);

% Calculate h and e for every integration step
for k = 1:l
    r = Yout_p (1:3, k);
    v = Yout_p (4:6, k);
    h_vect_p(:, k) = cross(r,v); 
    h_p(k) = norm( h_vect_p(:, k) );

    e_vect_p(:, k) = 1/mu_E*cross(v, h_vect_p(:, k))-r/norm(r);
    e_p(k) = norm(e_vect_p(:,k));
end

figure
plot(Tout_p, h_vect_p(1,:), 'b-.', Tout_p, h_vect_p(2,:), 'r-.', Tout_p, h_vect_p(3, :), 'g-.', Tout_p, h_p, 'k-.', 'LineWidth', 2);
legend('h_x', 'h_y', 'h_z', '|h|');
xlabel ('t [s]');
ylabel ('h_x, h_y, h_z, |h|  [km^2/s]');
title('Angular momentum');
grid on;

figure
plot(Tout_p, e_vect_p(1,:), 'b-.', Tout_p, e_vect_p(2,:), 'r-.', Tout_p, e_vect_p(3, :), 'g--', Tout_p, e_p, 'k--', 'LineWidth', 2);
legend('e_x', 'e_y', 'e_z', '|e|');
xlabel ('t [s]');
ylabel ('e_x, e_y, e_z, |e|  [-]');
title('Eccentricity');
grid on;

%% Perpendicularity check

% We want now to check that h and e remain perpendicular during each step
% Let's calculate scalar product and plot the error, it should be
% confrontabile con la tolleranza

err_perp = dot( h_vect_p, e_vect_p );

figure
plot(Tout, err_perp, 'b-')
legend('h dot e');
xlabel ('t [s]');
ylabel ('h dot e  [km^2/s]');
title('e-h dot product');
grid on;

%% Specific energy plot

% Dimension definition
eps_p = zeros(1, l);

for k = 1:l
    r = Yout_p (1:3, k);
    v = Yout_p (4:6, k);

    vnorm = norm(v);
    rnorm = norm(r);

    eps_p(k) = vnorm^2/2 - mu_E/rnorm;
end

figure
plot(Tout_p, eps_p, 'r-'),
legend('eps');
xlabel ('t [s]');
ylabel ('eps [km^2/s^2]');
title('Energy');
grid on;

%% Radial and transversal velocity plot

% Dimension definition
vr_p = zeros(1, l); % radial velocity
vt_p = zeros(1, l); % transversal velocity

for k = 1:l
    r = Yout_p (1:3, k);
    v = Yout_p (4:6, k);

    ur = r/norm(r); % radial unit vector
    vr_p(k) = dot(v, ur);

    uh = h_vect_p(:,k)/norm(h_vect_p(:,k)); % out-of-plane unit vector
    ut = cross(uh, ur); %transversal unit vector
    vt_p(k) = dot(v, ut);
end

figure
plot(Tout_p, vr_p, 'b-', Tout_p, vt_p, 'g-', 'LineWidth', 2)
legend('v_r', 'v_t');
xlabel ('t [s]');
ylabel ('v_r, v_t [km/s]');
title('Radial and transversal velocity');
grid on;

