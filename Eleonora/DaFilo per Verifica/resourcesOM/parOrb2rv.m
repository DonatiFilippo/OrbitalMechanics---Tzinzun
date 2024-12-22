function [rr, vv] = parOrb2rv(a, e, i, OM, om, theta, mu)

% Transformation from orbital elements to Cartesian state
%
% [rr, vv] = parorb2rv(a, e, i, OM, om, theta)
%
% Input arguments:
% --------------------------------------------------------
% a      [1x1]  semi-major axis     [km]
% e      [1x1]  eccentricity        [-]
% i      [1x1]  inclination         [rad]
% OM     [1x1]  RAAN                [rad]
% om     [1x1]  pericenter anomaly  [rad]
% theta  [1x1]  true anomaly        [rad]
% mu     [1x1]  gravitational parameter  [km^3/s^2]
% 
% --------------------------------------------------------
% Output arguments
% rr     [3x1]  position vector     [km]
% vv     [3x1]  velocity vector     [km/s]

% --------------------------------------------------------
%%
if nargin < 7
    mu = 398600;
end

% Rotation Matrix Definition:
R1 = [cos(OM), sin(OM), 0;
     -sin(OM), cos(OM), 0;
         0   ,   0    , 1];

R2 = [1,    0   ,   0   ;
      0,  cos(i), sin(i);
      0, -sin(i), cos(i)];

R3 = [cos(om), sin(om), 0;
     -sin(om), cos(om), 0;
         0   ,    0   , 1];

Tpf = (R3 * R2 * R1)';

% Orbital Parameters Evaluation:
p = a * (1 - (e^2));
h = sqrt(p * mu);

% Velocity and position vectors evaluation:
rrpf = [((h^2) / mu * 1 / (1+e * cos(theta))) * cos(theta);  
        ((h^2) / mu * 1 / (1+e * cos(theta))) * sin(theta); 
                                    0                     ];

vvpf = [mu / h * (-sin(theta))   ; 
        mu / h * (e + cos(theta)); 
                    0            ];

rr = Tpf * rrpf;
vv = Tpf * vvpf;

end