function dy = ode_2bp_perturbed(t, y, mu, muM, R, J2, date)
%eq_motion ODE system for the perturbed two-body problem (Keplerian motion)
% takes into account the effect of the second zonal harmonic J2 and the
% gravitational influence of the Moon
%
% PROTOTYPE:
% dy = eq_motion(t, y, mu, R, J2)
%
% INPUT:
% y[6x1]    State of the body ( rx, ry, rz, vx, vy, vz)        [L, L/T]
% mu[1]     Gravitational parameter of the primary             [L^3/T2]
% aJ2[3x1]  Perturbing acceleration due to J2                  [L/T^2]
% a3B[3x1]  Perturbing acceleration due to teh Moon            [L/T^2]
%
% OUTPUT: 
% dy[6X1]  Derivative of the state [L/T, L/T^2]
%
% CONTRIBUTORS:
%  Eleonora Domenichelli
%
% VERSION:
%  2024/10/14: First version
%  2024/12/22: Second version (addition of perturbing accelerations and
%  structural changes of the function)
%-------------------------------------------------------------------------

% Position
r = y(1:3);

% Distance from the primary
r_norm = norm(r);

% Derivative of the state
dy = [ y(4:6)
      (-mu/r_norm^3)*r + aJ2(r, mu, R, J2) + a3B_M(t, r, muM, date)];
end

% The time variable is needed to perform the integration?