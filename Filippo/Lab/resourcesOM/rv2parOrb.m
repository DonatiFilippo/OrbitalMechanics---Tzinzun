function [a, e, i, OM, om, theta] = rv2parOrb(rr, vv, mu)

% Transformation from Cartesian state to orbital elements
%
% [a, e, i, OM, om, theta] = rv2parorb(r, v, mu)
%
% Input arguments:
% --------------------------------------------------------
% rr     [3x1]  position vector     [km]
% vv     [3x1]  velocity vector     [km/s]
% mu     [1x1]  gravitational parameter  [km^3/s^2]
% 
% --------------------------------------------------------
% Output arguments
% a      [1x1]  semi-major axis     [km]
% e      [1x1]  eccentricity        [-]
% i      [1x1]  inclination         [rad]
% OM     [1x1]  RAAN                [rad]
% om     [1x1]  pericenter anomaly  [rad]
% theta  [1x1]  true anomaly        [rad]
% --------------------------------------------------------

    if nargin == 2
        mu = 398600;
    end
    
    % Needed Costants:
    K = [0; 0; 1];
    r = norm(rr);
    v = norm(vv);
    
    % Semimajor axis (a) evaluation:
    a = - (mu * r) / ((r * (v^2)) - 2 * mu);
    
    % Plane incination (i) evaluation:
    hh = cross(rr, vv);
    i = acos(hh(3, 1) / norm(hh));
    
    % RAAN (OM) evaluation:
    N = cross(K,hh);
    
    if N(2, 1) >= 0
        OM = acos(N(1, 1) / norm(N));     
    else
        OM = 2 * pi - acos(N(1, 1) / norm(N));  
    end
    
    % Eccentricity (e) evaluation:
    vr = dot(rr, vv) / r;
    ee = (1 / mu) * (((v^2) - (mu / r)) * rr - (r * vr * vv));
    e = norm(ee);
    
    % Argument of periapsis (om) evaluation:
    if ee(3, 1) >= 0
        om = acos(dot(N, ee) / (norm(N) * e));
    else
        om = 2 * pi - acos(dot(N, ee) / (norm(N) * e));
    end

    % True anomaly (theta) evaluation:
    if vr >= 0
        theta = acos(dot(ee, rr) / (e * r));
    else
        theta = 2 * pi - acos(dot(ee, rr) / (e * r));
    end

end


