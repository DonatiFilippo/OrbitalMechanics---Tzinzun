%% Eccentric anomaly in a specific point
clear;
clc;
close all;

% Datas
a = 7000; % [km]
mu_E = astroConstants(13); % [km^3/s^2]
E0 = 0;
t0 = 0;
T = 2*pi*sqrt(a^3/mu_E);
N = 99;
t = linspace(0, 2*T, N);
e = 0;
E = EccentricAnomalyScalar(t(end), a, e, mu_E, t0, E0);

%% Eccentric anomaly for different orbits 
clear;
clc;
close all;

% Datas
a = 7000; % [km]
mu_E = astroConstants(13); % [km^3/s^2]
E0 = 0;
t0 = 0;
T = 2*pi*sqrt(a^3/mu_E);
k = 2;
N = 100;
t = linspace(0, k*T, N);
e = [0, 0.2, 0.4, 0.6, 0.8, 0.95];

% Compute only one time e dimension
l = length(e);

% Initialize eccentric anomaly matrix and calculate for different orbits
E = zeros(l, N);

for i=1:l
    E(i, :) = EccentricAnomalyVector(t, a, e(i), mu_E, t0, N, E0);
    E(i, :) = rad2deg(E(i, :));
end

figure
hold on;
plot(t, E(1,:), 'b-', t, E(2, :), 'r--', t, E(3, :), 'y.', t, E(4, :), 'm-.',...
    t, E(5, :), 'g-', t, E(6,:), 'k--','LineWidth', 2)
grid on;
xlabel('t [s]')
ylabel('E [deg]')
legend('e = 0', 'e = 0.2', 'e = 0.4', 'e = 0.6', 'e = 0.8', 'e = 0.95');

% try 3D plot
e = [zeros(1, N); 0.2*ones(1, N); 0.4*ones(1, N); 0.6*ones(1, N);... 
    0.8*ones(1, N); 0.95*ones(1, N)];
t = t.*ones(6, N);

figure
surf(t, e, E)
xlabel('t [s]')
ylabel('e [-]')
zlabel('E [deg]')

