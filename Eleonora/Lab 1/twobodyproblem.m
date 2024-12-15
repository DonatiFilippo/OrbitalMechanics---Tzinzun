% TWO-BODY PROBLEM
clear;
clc;
close all;

% Parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
Rt = astroConstants(23); % Earth's radius

% Initial conditions
% HIGHLY ECCENTRIC AND INCLINED ORBIT
r0 = [6495; -970; -3622]; % [km]
v0 = [4.752; 2.130; 7.950]; % [km/s]
y0 = [r0; v0]; % [6x1]

% QUASI-CIRCULAR MEDIUM EARTH ORBIT
% r0 = [26578.137; 0; 0];
% v0 = [0; 2.221; 3.173]; 

% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [s]
tspan = linspace( 0, 5*T, 1000);
l = length(tspan);

% Options for the ODE solver
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

% Perform the integration
[Tout, Yout] = ode113 (@(t,y) ode_2bp(t,y,mu_E), tspan, y0, options);
Tout = Tout'; % [1xl] row vector
Yout = Yout'; % [6xl] row vector

% Plot the results
Earth_3D(Rt)
hold on
plot3( Yout(1, :), Yout(2, :), Yout(3, :), 'LineWidth', 2)
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;

%% Angular momentum and eccentricity plots

% We want to check that h and e remain constant in magnitude and direction 

% Dimensions definitions
h_vect = zeros(3, l); 
h = zeros(1, l); 
e_vect = zeros(3, l);
e = zeros(1, l);

% Calculate h and e for every integration step
for k = 1:l
    r = Yout (1:3, k);
    v = Yout (4:6, k);
    h_vect(:, k) = cross(r,v); 
    h(k) = norm( h_vect(:, k) );

    e_vect(:, k) = 1/mu_E*cross(v, h_vect(:, k))-r/norm(r);
    e(k) = norm(e_vect(:,k));
end

figure
plot(Tout, h_vect(1,:), 'b--', Tout, h_vect(2,:), 'r--', Tout, h_vect(3, :), 'g--', Tout, h, 'k--', 'LineWidth', 2);
legend('h_x', 'h_y', 'h_z', '|h|');
xlabel ('t [s]');
ylabel ('h_x, h_y, h_z, |h|  [km^2/s]');
title('Angular momentum');
grid on;

figure
plot(Tout, e_vect(1,:), 'b-', Tout, e_vect(2,:), 'r-', Tout, e_vect(3, :), 'g--', Tout, e, 'k--', 'LineWidth', 2);
legend('e_x', 'e_y', 'e_z', '|e|');
xlabel ('t [s]');
ylabel ('e_x, e_y, e_z, |e|  [-]');
title('Eccentricity');
grid on;

%% Perpendicularity check

% We want now to check that h and e remain perpendicular during each step
% Let's calculate scalar product and plot the error, it should be
% confrontabile con la tolleranza

err_perp = dot( h_vect, e_vect );

figure
plot(Tout, err_perp, 'b-')
legend('h dot e');
xlabel ('t [s]');
ylabel ('h dot e  [km^2/s]');
title('e-h dot product');
grid on;

%% Specific energy plot

% Dimension definition
eps = zeros(1, l);

for k = 1:l
    r = Yout (1:3, k);
    v = Yout (4:6, k);

    vnorm = norm(v);
    rnorm = norm(r);

    eps(k) = vnorm^2/2 - mu_E/rnorm;
end

figure
plot(Tout, eps, 'r-'),
legend('eps');
xlabel ('t [s]');
ylabel ('eps [km^2/s^2]');
title('Energy');
grid on;

%% Radial and transversal velocity plot

% Dimension definition
vr = zeros(1, l); % radial velocity
vt = zeros(1, l); % transversal velocity

for k = 1:l
    r = Yout (1:3, k);
    v = Yout (4:6, k);

    ur = r/norm(r); % radial unit vector
    vr(k) = dot(v, ur);

    uh = h_vect(:,k)/norm(h_vect(:,k)); % out-of-plane unit vector
    ut = cross(uh, ur); %transversal unit vector
    vt(k) = dot(v, ut);
end

figure
plot(Tout, vr, 'b-', Tout, vt, 'g-', 'LineWidth', 2)
legend('v_r', 'v_t');
xlabel ('t [s]');
ylabel ('v_r, v_t [km/s]');
title('Radial and transversal velocity');
grid on;