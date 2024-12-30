function [a, e, i, OM, om, theta] = rv2parorb(rr, vv, mu)

% rv2parorb - Convert Cartesian state vectors to Keplerian elements.
%
% PROTOTYPE: 
%   [a, e, i, OM, om, theta] = rv2parorb(rr, vv, mu)
%
% DESCRIPTION:
%   The function converts position and velocity vectors in Cartesian
%   coordinates to Keplerian elements.
%
% INPUT:
%   rr [3x1]      Position vector in Cartesian coordinates    [km]
%   vv [3x1]      Velocity vector in Cartesian coordinates    [km/s] 
%   mu [1x1]      Gravitational parameter of the primary      [km^3/s^2]
%
% OUTPUT:
%   a [1x1]       Semi-major axis                             [km]
%   e [1x1]       Eccentricity                                [-]
%   i [1x1]       Inclination                                 [rad]
%   OM [1x1]      Right ascension of ascending node           [rad]
%   om [1x1]      Pericenter anomaly                          [rad]
%   theta [1x1]   True anomaly                                [rad]
%
% CONTRIBUTORS:
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo 
%   Domenichelli Eleonora
%
%-------------------------------------------------------------------------

% Unit vectors of Geocentric Equatorial reference frame
I=[1; 0; 0];
J=[0; 1; 0];
K=[0; 0; 1];

% Magnitude of position (rr) and velocity (vv) vectors
r=norm(rr);
v=norm(vv);

% Radial component of velocity
v_r=dot(rr,vv)/r;

% Semi-major axis obtained from energy definition
a=1/(2/r-(v^2)/mu);

% Specific angular momentum vector and magnitude computation
hh=cross(rr,vv);
h=norm(hh);

% Inclination computation
i=acos(dot(hh,K)/h);

% Node line vector and magnitudde computation, considering the particular
% case where inclination is 0 (so OM will be zero)
if (i == 0) 
    NN=I;
else
    NN=cross(K,hh);
end
N=norm(NN);

% Right ascension of ascending node
if(dot(NN,J)>=0)
    OM=acos(dot(NN,I)/N);
else
    OM=2*pi-acos(dot(NN,I)/N);
end

% Eccentricity vector and magnitude computation, considering the particular
% case of circular orbits (e = 0)
ee=(cross(vv,hh))/mu-(rr/r);

if (norm(ee) <= 1e-8)
    ee=NN;
end
e=norm(ee);

% Pericenter anomaly
if(dot(ee,K)>=0)
    om=acos(dot(NN,ee)/(N*e));
else
    om=2*pi-acos(dot(NN,ee)/(N*e));
end

% True anomaly
if(v_r>=0)
    theta=acos(dot(ee,rr)/(e*r));
else
    theta=2*pi-acos(dot(ee,rr)/(e*r));
end
   
end