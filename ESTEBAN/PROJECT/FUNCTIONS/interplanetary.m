function dv = interplanetary (t1, t2, t3)
%% FUNCTION DESCRIPTION
%
%--------------------------------------------------------------------------
% DESCRIPTION:
% This code provides dv for a mission with a Lambert arc into a fly-by
% into a second Lambert arc
%
%--------------------------------------------------------------------------
% INPUTS:
%   t1          [1x1]       Time of departure           [mjd2000]
%   t2          [1x1]       Time of fly-by              [mjd2000]
%   t3          [1x1]       Time of arrival             [mjd2000]
%
%--------------------------------------------------------------------------
% OUTPUT:
%   dv          [1x1]       Total cost                  [km/s]
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
[kep1, ~] = uplanet(t1, 1); 
[kep2, ~] = uplanet(t2, 3); 
[kep3, ~, ~] = ephAsteroids(t3, 30);

muS = astroConstants(4);
muE = astroConstants(13); 

[r1, v1i] = kep2cart(kep1(1), kep1(2), kep1(3), kep1(4), kep1(5), kep1(6), muS); % Departure planet
[r2, v2f] = kep2cart(kep2(1), kep2(2), kep2(3), kep2(4), kep2(5), kep2(6), muS); % Fly-by planet
[r3, v3f] = kep2cart(kep3(1), kep3(2), kep3(3), kep3(4), kep3(5), kep3(6), muS); % Asteroid

%% Set time quantities
dt1 = (t2 - t1)*86400;
dt2 = (t3 - t2)*86400;

%% Solver
% Transfer leg 1
[~, ~, ~, ERROR1, v1t, v2t, ~, ~] = lambertMR(r1, r2, dt1, muS, 0, 0, 0, 0);

% Transfer leg 2
[~, ~, ~, ERROR2, v2t_1, v3t, ~, ~] = lambertMR(r2, r3, dt2, muS, 0, 0, 0, 0);

if ERROR1 == 0 && ERROR2 == 0
        vinfmin_vec = v2t' - v2f;
        vinfplus_vec = v2t_1' - v2f;
        [~, ~, ~, rp_fb, ~, ~, ~, ~, ~, ~, ~, ~, dv_fb, ~] = flyby_powered(vinfmin_vec, vinfplus_vec, muE);

        if not(isnan(rp_fb))
            dv_dep = norm(v1t' - v1i);
            dv_arr = norm(v3t' - v3f);
            dv = dv_dep + dv_fb + dv_arr;
        else
            dv = 2000; % Supposed to be equal to NaN, put equal to 2000 to not take it in account in minimization
        end
end

end
