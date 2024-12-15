function [E] = EccentricAnomalyScalar(tspan, a, e, mu, t0, E0)
%EccentricAnomalyScalar provides the value of eccentric anomaly E for different t
%in case of elliptic orbits only (0 <= e < 1) and with t0, E0 generic (!=0)
%
% PROTOTYPE:
% E = EccentricAnomalyScalar(t, a, e, mu, t0, E0)
%
% INPUT
% t[1]    Time                                              [T]
% a[1]    Semi-major axis of the orbit                      [L]
% e[1]    Eccentricity magnitude                            [-]
% mu[1]   Gravitational parameter of the primary            [L^3/T^2]
% t0[1]   Reference initial time                            [T]
% E0[1]   Reference initial eccentric anomaly               [rad] 
%
% OUTPUT:
% E[1]    Eccentric anomaly (function of time)              [rad]  
%
% CONTRIBUTORS:
%  Eleonora Domenichelli
%
% VERSION:
%  2024/10/15: First version
%
%-------------------------------------------------------------------------

% Define n
n = sqrt(mu/(a^3));

% Initial condition for eccentric anomaly considered optional
if nargin < 6
    E0 = n*t0 + (e*sin(n*t0) / (1 - sin(n*t0 + e) + sin(n*t0)));
end

% Define orbit's period
T = 2*pi/n;

k = floor((tspan-t0)/T);

% Find time remainder
% needed only if t is bigger than orbit's period, so we distinguish those
% two cases
if tspan-t0 > T
   rem = (tspan-t0)-k*T;
else
   rem = tspan-t0;
end

% Define Kepler's equation as anonymus function
kep = @(E) -E+e*sin(E)+n*rem+E0-e*sin(E0);

% Solve Kepler's equation
options = optimoptions('fsolve', 'Display','off');
E_rem = fsolve( kep , E0, options);


E = k*2*pi + E_rem;

end