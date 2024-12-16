mu = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
J = astroConstants(9); % Second Armonic Constant [-] 
R = astroConstants(23); % Earth's Radius [km]

% Initial Conditions
rr0 = [7571; 0.01; deg2rad(87.9)]; % Initial Position [Km]
vv0 = [deg2rad(180); deg2rad(180); deg2rad(0) ]; % Initial Velocity [Km/s^2]
s0 = vertcat(rr0, vv0);

%% Main Section

% Preliminary Calculations
a = s0(1); % Semi-major axis [km]
 
% Settings
T = 2*pi*sqrt( a^3/mu ); % Orbital period [s]
tspan = linspace( 0, 100*T, 100000 );
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

% Integrtor
[ts, s] = ode113(@(t, y) gaussFun(t, y, acc_pert_fun( y, J, mu, R), mu), tspan, s0, options);


figure();
hold on
plot(ts/T,s(:,1));

%% Filtering

a = s(:, 1);
a_fil = movmean(a, 1000, 'Endpoints','fill');

plot(ts/T, a_fil);
