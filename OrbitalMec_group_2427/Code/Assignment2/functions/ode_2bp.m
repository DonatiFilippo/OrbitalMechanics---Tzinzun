function dy = ode_2bp(~, y, mu_E)

% ode_2bp - ODE system for the two-body problem
%
% PROTOTYPE:
%  dy = ode_2bp(t, y, mu)
%
% DESCRIPTION:
%   The function gives the ODE system for the not perturbed two-body
%   problem.
%
% INPUT:
% t[1]     Time (can be omitted as the system is autonomous)  [s]
% y[6x1]   State of the body ( rx, ry, rz, vx, vy, vz)        [km, km/s]
% mu[1]    Gravitational parameter of the primary             [km^3/s^2]
%
% OUTPUT: 
% dy[6X1]  Derivative of the state                            [km/s, km/s^2]
%
% CONTRIBUTORS:
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo 
%   Domenichelli Eleonora
%
%-------------------------------------------------------------------------

% Position
r_vect= y(1:3);

% Distance from the primary
r=norm(r_vect);

% Derivative of the state
dy = [ y(4:6)
      (-mu_E/r^3)*y(1:3)];
end