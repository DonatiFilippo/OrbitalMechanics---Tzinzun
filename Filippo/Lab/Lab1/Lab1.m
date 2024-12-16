%% Constant And Initial Conditions

% Constants
mu = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
J = 0.00108263; % Second Armonic Constant [-] 
R = astroConstants(23); % Earth's Radius [km]

% Initial Conditions
rr0 = [26578.137; 0; 0]; % Initial Position [Km]
vv0 = [0; 2.221; 3.173 ]; % Initial Velocity [Km/s^2]
s0 = vertcat(rr0, vv0);

%% Main Section

% Preliminary Calculations
a = 1/( 2/norm(rr0) - dot(vv0,vv0)/mu ); % Semi-major axis [km]
 
% Settings
T = 2*pi*sqrt( a^3/mu ); % Orbital period [s]
tspan = linspace( 0, 2*T, 1000 );
options2 = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Integrtor
[ts, s] = ode113(@(t, y) twoBodyProblemPert(t, y, mu, J, R), tspan, s0, options2);

%% Plot

% Orbital Plot
figure()
plot3( s(:,1), s(:,2), s(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;