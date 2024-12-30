function A_RSW = TNH2RSW(kep, mu)

% TNH2RSW - Rotating matrix from TNH to RSW reference frame
%
% PROTOTYPE: 
%   A_RSW = TNH2RSW(kep, mu)
%
% DESCRIPTION:
%   The function gives the rotation matrix from TNH (tangential-normal-out
%   of plane) reference framae to RSW (radial-transversal-out of plane)
%   reference frame, taking as input the Keplerian elements vector.
%   The first frame {t, n, h} is characterised by:
%       t-axis: on the orbit plane, along the direction of the velocity
%       h-axis: out of plane, along the direction of the angular momentum
%       of the orbit
%       n-axis: on the orbit plane, completes the reference frame
%   The second frame {r, w, s} is characterized by:
%       r-axis: on the orbit plane, along the direction of the radius
%       w-axis: out of plane, along the direction of the angular momentum
%       of the orbit
%       s-axis: on the orbit plane, completes the reference frame
%
% INPUT:
%   kep [1x6]    Keplerian parameters vector
%                (a, e, i, OM, om, theta)                  [km, -, rad]
%   mu [1x1]     Gravitational parameter of the primary    [km^3/s^2]
%
% OUTPUT: 
%   A_RSW [3x3]   Rotation matrix                          [-]
%
% CONTRIBUTORS:
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo 
%   Domenichelli Eleonora
%
%-------------------------------------------------------------------------

% Parameters extraction
a = kep(1);
e = kep(2);
theta = kep(6);

% Important terms computation
p = a*(1-e^2);
r = p / ( 1+e*cos(theta));
v = sqrt(2*mu/r - mu/a);
n = sqrt(mu/(a^3));
b = a*sqrt(1-e^2);
h = n * a * b;

% Rotation matrix
A_RSW = h/(p*v)*[e*sin(theta), -(1+e*cos(theta)), 0;
                  1+e*cos(theta), e*sin(theta), 0;
                  0, 0, (p*v)/h];
end