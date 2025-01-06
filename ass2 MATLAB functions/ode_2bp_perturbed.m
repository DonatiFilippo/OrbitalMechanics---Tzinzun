function dy = ode_2bp_perturbed(t, y, mu, muM, R, J2, date)

% ode_2bp_perturbed - ODE system for the perturbed two-body problem
%
% PROTOTYPE:
%   dy = eq_motion(t, y, mu, R, J2)
%
% DESCRIPTION:
%   The functions gives the ODE system for the perturbed two-body problem,
%   taking into account the effect of the second zonal harmonic J2 and the
%   gravitational influence of the Moon.
%
% INPUT:
%   t [1x1]    Time instant for acceleration evaluation          [s]
%   y [6x1]    State of the body (rx, ry, rz, vx, vy, vz)        [km, km/s]
%   mu [1x1]   Gravitational parameter of the primary            [km^3/s^2]
%   muM [1x1]  Gravitational parameter of the Moon               [km^3/s^2]
%   R [1x1]    Radius of the primary                             [km]
%   J2 [1x1]   Gravitational harmonic coefficient of the primary [-]
%
%   date [1x6]    Date in the Gregorian calendar, as a 6-elements vector
%                 [year, month, day, hour, minute, second]. For dates 
%                 before 1582, the resulting date components are valid 
%                 only in the Gregorian proleptic calendar. This is based
%                 on the Gregorian calendar but extended to cover dates 
%                 before its introduction. Date must be after 12:00 noon,
%                 24 November -4713.
%
% OUTPUT: 
%   dy[6X1]  Derivative of the state                              [km/s, km/s^2]
%
% CONTRIBUTORS:
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo 
%   Domenichelli Eleonora
%
%-------------------------------------------------------------------------

% Position extraction
r = y(1:3);

% Distance from the primary
r_norm = norm(r);

% Derivative of the state
dy = [ y(4:6)
      (-mu/r_norm^3)*r + aJ2(r, mu, R, J2) + a3B_M(t, r, muM, date)];
end