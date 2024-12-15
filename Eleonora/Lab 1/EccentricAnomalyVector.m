function [E_vect] = EccentricAnomalyVector(tspan, a, e, mu, t0, N, E0)
%EccentricAnomalyVector provides the value of eccentric anomaly E for different t
%in case of elliptic orbits only (0 <= e < 1)
%
% PROTOTYPE:
% E = EccentricAnomalyVector(t, a, e, mu, t0, N, k, E0)
%
% INPUT
% t[1xN]    Time                                              [T]
% a[1]      Semi-major axis of the orbit                      [L]
% e[1]      Eccentricity magnitude                            [-]
% mu[1]     Gravitational parameter of the primary            [L^3/T^2]
% t0[1]     Reference initial time                            [T]
% N[1]      Time vector length                                [-]
% E0[1]     Reference initial eccentric anomaly               [rad]  
%
% OUTPUT:
% E_vect[1xN]    Eccentric anomaly (function of time)              [rad]  
%
% CONTRIBUTORS:
%  Eleonora Domenichelli
%
% VERSION:
%  2024/10/17: First version
%
%-------------------------------------------------------------------------
% Define n
n = sqrt(mu/(a^3));

% Initial condition for eccentric anomaly considered optional
if nargin < 8
    E0 = n*t0 + (e*sin(n*t0) / (1 - sin(n*t0 + e) + sin(n*t0)));
end

% Define orbit's period
T = 2*pi/n;

% Initialize eccentric anomaly vector
E_vect = zeros(1,N);

for j = 1:N

k = floor((tspan(j)-t0)/T);
% Find time remainder
% needed only if t is bigger than orbit's period, so we distinguish those
% two cases
if tspan(j)-t0 > T
   rem = (tspan(j))-k*T;
else
   rem = tspan(j)-t0;
end

% Define Kepler's equation as anonymus function
kep = @(E) -E+e*sin(E)+n*rem+E0-e*sin(E0);

% Solve Kepler's equation
options = optimoptions('fsolve', 'Display','off');
E_rem = fsolve( kep , E0, options);


E_vect(1,j) = k*2*pi + E_rem;
end

end