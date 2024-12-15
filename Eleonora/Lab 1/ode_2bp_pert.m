function dy = ode_2bp_pert(~, y, mu_E, Rt, J2)
%ode_2bp ODE system for the perturbed two-body problem (Keplerian motion)
%thakes into account the effect of the second zonal harmonic J2
%
% PROTOTYPE:
%  dy = ode_2bp(t, y, mu, R, J2)
%
% INPUT:
% t[1]     Time (can be omitted as the system is autonomous)  [T]
% y[6x1]   State of the body ( rx, ry, rz, vx, vy, vz)        [L, L/T]
% mu[1]    Gravitational parameter of the primary             [L^3/T2]
% R[1]     Primary's radius                                   [L]
% J2[1]    Second zonal harmonic                              [-]
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

% Perturbation
a = 3/2*( J2*mu_E*Rt^2 )/r^4 * [r_vect(1)/r * (5*(r_vect(3)/r)^2-1);...
                                r_vect(2)/r * (5*(r_vect(3)/r)^2-1);...
                                r_vect(3)/r * (5*(r_vect(3)/r)^2-3)];

% Derivative of the state
dy = [ y(4:6)
      (-mu_E/r^3)*y(1:3) + a];
end