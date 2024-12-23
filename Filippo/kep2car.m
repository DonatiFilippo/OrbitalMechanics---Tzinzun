function [rr,vv] = kep2car(a, e, i, OM, w, f, mu)
%KEP2CAR Summary of this function goes here
%   Detailed explanation goes here

% Rotation Matrix Definition:
R1 = [cos(OM), sin(OM), 0;
     -sin(OM), cos(OM), 0;
         0   ,   0    , 1];

R2 = [1,    0   ,   0   ;
      0,  cos(i), sin(i);
      0, -sin(i), cos(i)];

R3 = [cos(w), sin(w), 0;
     -sin(w), cos(w), 0;
         0   ,    0   , 1];

Tpf = (R3 * R2 * R1)';

p = a * (1 - (e^2));
h = sqrt(p * mu);

% Velocity and position vectors evaluation:
rrpf = [((h^2) / mu * 1 / (1+e * cos(f))) * cos(f);  
        ((h^2) / mu * 1 / (1+e * cos(f))) * sin(f); 
                                    0              ];

vvpf = [mu / h * (-sin(f))   ; 
        mu / h * (e + cos(f)); 
                    0            ];

rr = Tpf * rrpf;
vv = Tpf * vvpf;

end

