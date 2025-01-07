function [dv, dv_dep, dv_arr, r1, v1i, r2, v2f, r3, v3f, v1t, v2t, v2t_1, v3t, vinfmin_vec, vinfplus_vec, cost] = interplanetary (t1, t2, t3)
%% FUNCTION DESCRIPTION
%
%--------------------------------------------------------------------------
% DESCRIPTION:
% This code provides dv for a mission with a Lambert arc into a fly-by
% into a second Lambert arc
%
%--------------------------------------------------------------------------
% INPUTS:
%
%   t1          [1x1]       Time of departure           [mjd2000]
%   t2          [1x1]       Time of fly-by              [mjd2000]
%   t3          [1x1]       Time of arrival             [mjd2000]
%
%--------------------------------------------------------------------------
% OUTPUT:
%   dv          [1x1]       Total cost                  [km/s]
%   dv_dep      [1x1]       Delta-v for departure arc   [km/s]
%   dv_arr      [1x1]       Delta-v for arrival arc     [km/s]
%   r1          [3x1]       Initial position vector     [km]
%   v1i         [3x1]       Initial velocity vector     [km/s]
%   r2          [3x1]       Fly-by position vector      [km]
%   v2f         [3x1]       Fly-by velocity vector      [km/s]
%   r3          [3x1]       Arrival position vector     [km]
%   v3f         [3x1]       Arrival velocity vector     [km/s]
%   v1t         [3x1]       Transfer velocity vector at departure [km/s]
%   v2t         [3x1]       Transfer velocity vector at fly-by (arrival) [km/s]
%   v2t_1       [3x1]       Transfer velocity vector at fly-by (departure) [km/s]
%   v3t         [3x1]       Transfer velocity vector at arrival [km/s]
%   vinfmin_vec [3x1]       Velocity difference at infinity (approaching fly-by planet) [km/s]
%   vinfplus_vec[3x1]       Velocity difference at infinity (leaving fly-by planet) [km/s]
%   cost        [1x1]       Penalty cost for infeasible solutions (large value or NaN)
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

%% Initial condition

if t2 <= t1 || t3 <= t2
        cost = 1e6; 
         dv = NaN;
else
    cost = NaN;

[kep1, ~] = uplanet(t1, 1); 
[kep2, ~] = uplanet(t2, 3); 
[kep3, ~, ~] = ephAsteroids(t3, 30);

muS = astroConstants(4);
muE = astroConstants(13); 

[r1, v1i] = kep2cart(kep1(1), kep1(2), kep1(3), kep1(4), kep1(5), kep1(6), muS); % Departure planet
[r2, v2f] = kep2cart(kep2(1), kep2(2), kep2(3), kep2(4), kep2(5), kep2(6), muS); % Fly-by planet
[r3, v3f] = kep2cart(kep3(1), kep3(2), kep3(3), kep3(4), kep3(5), kep3(6), muS); % Asteroid

dt1 = (t2 - t1)*86400;
dt2 = (t3 - t2)*86400;

% Transfer leg 1
[~, ~, ~, ERROR1, v1t, v2t, ~, ~] = lambertMR(r1, r2, dt1, muS, 0, 0, 0, 0);

% Transfer leg 2
[~, ~, ~, ERROR2, v2t_1, v3t, ~, ~] = lambertMR(r2, r3, dt2, muS, 0, 0, 0, 0);

if ERROR1 == 0 && ERROR2 == 0
        vinfmin_vec = v2t' - v2f;
        vinfplus_vec = v2t_1' - v2f;
        [~, ~, ~, rp_fb, ~, ~, ~, ~, ~, ~, ~, ~, ~, dv_fb] = flyby_powered(vinfmin_vec, vinfplus_vec, muE);

        if not(isnan(rp_fb))
            dv_dep = norm(v1t' - v1i);
            dv_arr = norm(v3t' - v3f);
            dv = dv_dep + dv_fb + dv_arr;
        else
            dv = NaN;
        end
end
end
end
