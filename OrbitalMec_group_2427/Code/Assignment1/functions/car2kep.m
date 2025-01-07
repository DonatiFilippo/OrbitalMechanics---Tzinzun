function [a, e, i, Omega, omega, nu] = car2kep(varargin)
%% FUNCTION DESCRIPTION
%
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function converts a state vector (position and velocity in Cartesian 
% coordinates) into Keplerian orbital elements. It is applicable for any orbit 
% around a central body, assuming a Newtonian two-body problem.
%
%--------------------------------------------------------------------------
% INPUTS:
%   varargin    [1x6] or [1x2]
%       - A single input vector containing position and velocity components:
%         [r_x, r_y, r_z, v_x, v_y, v_z] [km, km/s].
%       - Alternatively, the position and velocity vectors can be provided 
%         separately as [r] and [v] ([km] and [km/s]).
%
%--------------------------------------------------------------------------
% OUTPUTS:
%   a       [1x1]   Semi-major axis                                [km]
%   e       [1x1]   Eccentricity                                   [-]
%   i       [1x1]   Inclination                                    [deg]
%   Omega   [1x1]   Right ascension of the ascending node (RAAN)   [deg]
%   omega   [1x1]   Argument of periapsis                          [deg]
%   nu      [1x1]   True anomaly                                   [deg]
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

if nargin == 0 % if no inputs, break
    error('At least 1 input argument required.');
else    
    if nargin < 3 % if mu is not specified, use Earth's
        constants; % load constants from file
        mu = mu_E;
    end
    if nargin == 1 % if X is a state-vector and V is not specified
        if length(X) == 6
            V = X(4:6);
            X = X(1:3);
        else
            error('The state-vector has to have 6 components.');
        end
    else % check X and V have 3 components each
        if length(X) ~= 3 || length(V) ~= 3
            error('The position and velocity vectors have to have 3 components each.');
        end
    end
    
    r = norm(X);
    v = norm(V);
    h = cross(X,V);
    N = cross([0;0;1],h);
    
    a = 1/(2/r - v^2/mu);
    e = 1/mu*cross(V,h) - X/r;
    i = acos(h(3)/norm(h));
    Nxy = sqrt(N(1)^2 + N(2)^2);
    raan = atan2(N(2)/Nxy,N(1)/Nxy);
    
    NN = N/norm(N);
    ee = e/norm(e);
    omega = sign(dot(cross(NN,e),h))*acos(dot(ee,NN));
    f = sign(dot(cross(e,X),h))*acos(dot(X/r,ee));
    
    e = norm(e);
    i = i*180/pi;
    if isnan(raan) % If the raan is not defined, assign (arbitrary) 0 value
        raan = 0;
    else
        raan = raan*180/pi;
        if raan < 0
            raan = raan + 360;
        end
    end
    if isnan(omega) % If omega is not defined, assign (arbitrary) 0 value
        omega = 0;
    else
        omega = omega*180/pi;
        if omega < 0
            omega = omega + 360;
        end
    end
    f = f*180/pi;
    if f < 0
        f = f + 360;
    end
    M = meanAnomaly(f,e);
end
end

function M = meanAnomaly(f,e)
% function M = meanAnomaly(f,e)
% Computes mean anomaly for an elliptical orbit
% Required inputs:
%   f: true anomaly [deg]
%   e: eccentricity [-]
% Outputs:
%   M: mean anomaly [deg]

f = f*pi/180;
E = 2*atan(sqrt((1 - e)/(1 + e))*tan(f/2)); % compute eccentric anomaly
M = E - e*sin(E); % compute mean anomaly
M = M*180/pi;
if M < 0
    M = M + 360;
end
end
