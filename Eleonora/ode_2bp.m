function dy = ode_2bp(~, y, mu_E)
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE:
%  dy = ode_2bp(t, y, mu)
%
% INPUT:
% t[1]     Time (can be omitted as the system is autonomous)  [T]
% y[6x1]   State of the body ( rx, ry, rz, vx, vy, vz)        [L, L/T]
% mu[1]    Gravitational parameter of the primary             [L^3/T2]
%
% OUTPUT: 
% dy[6X1]  Derivative of the state [L/T, L/T^2]
%
% CONTRIBUTORS:
%  Eleonora Domenichelli
%
% VERSION:
%  2024/10/14: First version
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