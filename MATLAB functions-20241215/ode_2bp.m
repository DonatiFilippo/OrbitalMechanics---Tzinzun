function dy = ode_2bp( ~ , y , mu )
%% FUNCTION DESCRIPTION
%
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function defines the system of ordinary differential equations (ODEs) 
% for the two-body problem (Keplerian motion). It calculates the derivatives 
% of the state vector (position and velocity) under the influence of a central 
% gravitational force.
%
%--------------------------------------------------------------------------
% INPUTS:
%   ~         [1x1]   Time variable (not used since the system is autonomous) [-]
%   y         [6x1]   State vector (position and velocity):
%                     [rx, ry, rz, vx, vy, vz] [km, km/s]
%   mu        [1x1]   Gravitational parameter of the central body     [km^3/s^2]
%
%--------------------------------------------------------------------------
% OUTPUTS:
%   dy        [6x1]   Derivative of the state vector:
%                     [vx, vy, vz, ax, ay, az] [km/s, km/s^2]
%
%--------------------------------------------------------------------------
% Group number : 27
%
% Created and maintained by : 
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo 
%   Domenichelli Eleonora
%
%--------------------------------------------------------------------------

r = y(1 : 3);
v = y (4 : 6);

rnorm = norm(r);

dy = [v 
        (-mu/rnorm^3)*r];

end
