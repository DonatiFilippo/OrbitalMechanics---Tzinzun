function [vinfm, vinfp, delta, rp, am, ap, em, ep, vpm, vpp, deltam, deltap, dv, dvp] = flyby_powered(vinfmin_vec, vinfplus_vec, mu)
%% FUNCTION DESCRIPTION
%
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function computes the parameters and velocity changes associated with 
% a powered gravity-assist maneuver (flyby). It calculates the trajectory 
% characteristics and required delta-v for modifying the spacecraft's trajectory 
% during the flyby of a celestial body.
%
%--------------------------------------------------------------------------
% INPUTS:
%   vinfmin_vec  [1x3]   Incoming hyperbolic excess velocity vector before flyby [km/s]
%   vinfplus_vec [1x3]   Outgoing hyperbolic excess velocity vector after flyby  [km/s]
%   mu           [1x1]   Gravitational parameter of the flyby planet             [km^3/s^2]
%
%--------------------------------------------------------------------------
% OUTPUTS:
%   vinfm        [1x1]   Magnitude of incoming velocity at infinity              [km/s]
%   vinfp        [1x1]   Magnitude of outgoing velocity at infinity              [km/s]
%   delta        [1x1]   Turning angle during the flyby                          [rad]
%   rp           [1x1]   Pericenter radius (closest approach to the planet)      [km]
%   am           [1x1]   Semi-major axis of the incoming hyperbolic trajectory   [km]
%   ap           [1x1]   Semi-major axis of the outgoing hyperbolic trajectory   [km]
%   em           [1x1]   Eccentricity of the incoming hyperbolic trajectory      [-]
%   ep           [1x1]   Eccentricity of the outgoing hyperbolic trajectory      [-]
%   vpm          [1x1]   Pericenter velocity for the incoming trajectory         [km/s]
%   vpp          [1x1]   Pericenter velocity for the outgoing trajectory         [km/s]
%   deltam       [1x1]   Incoming hyperbolic trajectory deflection angle         [rad]
%   deltap       [1x1]   Outgoing hyperbolic trajectory deflection angle         [rad]
%   dv           [1x1]   Total delta-v required for the maneuver                 [km/s]
%   dvp          [1x1]   Delta-v applied at the pericenter                       [km/s]
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
%

vinfm = norm(vinfmin_vec);
vinfp = norm(vinfplus_vec);

delta = acos((dot(vinfmin_vec, vinfplus_vec))/(vinfm*vinfp));

rp_test = 6378.136; % Initial guess
options = optimset('Display','off');

fun = @(rp) -delta + asin(1/(1 + rp/mu*vinfm^2)) + asin(1/(1 + rp/mu*vinfp^2));

if ~isfinite(fun(rp_test))
    error('fzero value is not valid');
end

% Solver
try
    rp = fzero(fun, rp_test, options);
catch
    rp = NaN; 
end

h_condition = 800;
rp_critical = 6378.136 + h_condition;

if rp > rp_critical
    am = -mu/vinfm^2;
    ap = -mu/vinfp^2;
    em = 1+rp*vinfm^2/mu;
    ep = 1+rp*vinfp^2/mu;
    vpm = sqrt(2*mu*(1./rp + 1/2/abs(am)));
    vpp = sqrt(2*mu*(1./rp + 1/2/abs(ap)));
    deltam = 2*asin(1/em);
    deltap = 2*asin(1/ep);
    dv = norm(vinfplus_vec - vinfmin_vec);
    dvp = abs(vpp - vpm);
else
    rp = NaN;
    am = NaN;
    ap = NaN;
    em = NaN;
    ep = NaN;
    vpm = NaN;
    vpp = NaN;
    deltam = NaN;
    deltap = NaN;
    dv = NaN;
    dvp =NaN;
end

end
