% Orbit disturbances and propagation
clear
clc
close all

 % Set your initial time and state (Cartesian or Keplerian): t0, s0
 % Set your integration time span: tspan
 % Set physical parameters (e.g., gravitational parameter of the primary, J2, etc.)
 % Set ODE solver options (e.g., RelTol, AbsTol): options
 % Numerical integration of the equations of motion

% Initial conditions
i = 87.9*(pi/180);
a0 = 7571;
kep0 = [a0;0.01;i;pi;pi;0];
% Phisical parameters
mu_E = astroConstants(13);
Re = astroConstants(23);
J2 = astroConstants(9);
T = 2.*pi.*sqrt(a0^3/mu_E);
tmax = 100.*T;
tspan = linspace(0,tmax,1000);

% Integration options 
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Solver run
[Tsim1,KEP1] = ode113(@(t, kep) eq_motion(t, kep, acc_pert_fun(kep, J2, mu_E, Re), J2, mu_E, Re), tspan, kep0, options);

for k = 1:length(tspan)
    KEP1(k,3:6) = [rad2deg(KEP1(k,3)), rad2deg(KEP1(k,4)), rad2deg(KEP1(k,5)), rad2deg(KEP1(k,6))] ;
end


%% Cartesian coordinates propagator
% Integration of perturbed 2 body problem, considering earth oblateness J2


% Initial conditions
[r0,v0] = parOrb2rv(kep0(1),kep0(2),kep0(3),kep0(4),kep0(5),kep0(6));
y0 = [r0; v0];

% Integration options and calculations
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[Tsim2,Y] = ode113(@(t,y) ode_2bp_pertJ2(t,y,mu_E,Re,J2),tspan,y0,options);

KEP2 = zeros(length(tspan),6);
for k = 1:length(tspan)
    [a, e, i, OM, om, theta] = car2kep(Y(k,1:3),Y(k,4:6), mu_E);
    KEP2(k,1:6) = [a, e, rad2deg(i), rad2deg(OM), rad2deg(om), (theta)];
end
KEP2(:,6) = unwrap(KEP2(:,6));
KEP2(:,6) = rad2deg(KEP2(:,6));

%% Postprocess and results comparison
rel_err = [];
for k = 1:length(tspan)
    err = abs(KEP2(k,1)-KEP1(k,1))/a0;
    rel_err = [rel_err err];
end
Tsim1 = Tsim1./T;
semilogy(Tsim1,rel_err)
grid on
%% 
figure
Tsim1 = Tsim1/T;
plot(Tsim1,KEP1(:,1))
legend('Gauss Equations','Cartesian')
grid on