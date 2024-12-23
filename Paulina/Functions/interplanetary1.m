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
best_dv = inf; % Initialize the best delta-v
Nrev_max = 3; % Maximum number of revolutions to consider
for Nrev = 0:Nrev_max
    % Lambert arc from departure to flyby
    [~, ~, ~, ERROR1, v1t, v2t, ~, ~] = lambertMR(r1, r2, dt1, muS, 0, Nrev, 0, 0);
    % Lambert arc from flyby to arrival
    [~, ~, ~, ERROR2, v2t_1, v3t, ~, ~] = lambertMR(r2, r3, dt2, muS, 0, Nrev, 0, 0);

    % Check for valid solutions
    if ERROR1 == 0 && ERROR2 == 0
        % Compute flyby parameters
        vinfmin_vec = v2t' - v2f; % Incoming hyperbolic excess velocity
        vinfplus_vec = v2t_1' - v2f; % Outgoing hyperbolic excess velocity
        [~, ~, ~, rp_fb, ~, ~, ~, ~, ~, ~, ~, ~, dv_fb, ~] = flyby_powered(vinfmin_vec, vinfplus_vec, muE);

        % Check if flyby radius is feasible
        if ~isnan(rp_fb)
            % Compute delta-v contributions
            dv_dep = norm(v1t' - v1i); % Departure delta-v
            dv_arr = norm(v3t' - v3f); % Arrival delta-v
            dv_total = dv_dep + dv_fb + dv_arr;

            % Update the best delta-v
            if dv_total < best_dv
                best_dv = dv_total;
            end
        end
    end
end

%% Output the best delta-v
if best_dv < inf
    dv = best_dv; % Return the best delta-v
else
    dv = NaN; % No valid solution found
end

end