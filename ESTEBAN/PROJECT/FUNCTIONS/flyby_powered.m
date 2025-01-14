function [vinfm, vinfp, delta, rp, am, ap, em, ep, vpm, vpp, deltam, deltap, dv, dvp] = flyby_powered(vinfmin_vec, vinfplus_vec, mu)
%% FUNCTION DESCRIPTION
%
%--------------------------------------------------------------------------
% DESCRIPTION:
% This code provides all data correlated to a powered gravity-assist 
% manoeuvre
%
%--------------------------------------------------------------------------
% INPUTS:
%
% vinfmin_vec     [1x3]    Velocity at infinite before fb    [km/s]
% vinfplus_vec    [1x3]    Velocity at infinite after fb     [km/s]
%
%--------------------------------------------------------------------------
% OUTPUT:
%
%   TO DO ...
%
%--------------------------------------------------------------------------
% Group number : 27
%
% Created and maintained by : 
%
% Azevedo Da Silva Esteban
% Gavidia Pantoja Maria Paulina
% Donati Filippo 
% Domenichelli Eleonora
% 
%--------------------------------------------------------------------------
% LAST UPDATE: 21-12-2024
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