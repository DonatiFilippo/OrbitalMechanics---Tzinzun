function dvtot = Lambert (t1, t2, planet1, planet2)
%% FUNCTION DESCRIPTION
%
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function provides a fast and optimized implementation of the Lambert problem
% solver for interplanetary missions. It calculates the total delta-v required for 
% a transfer between two celestial bodies given their departure and arrival dates.
% The function utilizes the `LambertMR` solver, simplifying its application by 
% automatically retrieving planetary or asteroid ephemerides for the given dates.
%
%--------------------------------------------------------------------------
% INPUTS:
%   t1          [1x1]   Departure time (MJD2000)                      [mjd2000]
%   t2          [1x1]   Arrival time (MJD2000)                        [mjd2000]
%   planet1     [1x1]   ID of the departure celestial body            [-]
%                       (1 = Mercury, ..., 9 = Pluto, >10 = Asteroid ID)
%   planet2     [1x1]   ID of the arrival celestial body              [-]
%                       (1 = Mercury, ..., 9 = Pluto, >10 = Asteroid ID)
%
%--------------------------------------------------------------------------
% OUTPUT:
%   dvtot       [1x1]   Total delta-v required for the transfer [km/s]
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

 muS = astroConstants(4);

% Initial condition

if planet1 <= 10
    [kep1,~] = uplanet(t1, planet1);
else
    [kep1,~] = ephAsteroids(t1, planet1);
end

if planet2 <= 10
    [kep2,~] = uplanet(t2, planet2);
else
    [kep2,~] = ephAsteroids(t2, planet2);
end

muS = astroConstants(4);

[r1, v1i] = kep2cart(kep1, muS);
[r2, v2f] = kep2cart(kep2,muS);

% Set time quantities
dt = (t2 - t1)*24*60*60;

% Solver
[~, ~, ~, ~, v1t, v2t, ~, ~] = lambertMR(r1, r2, dt, muS, 0, 0, 0, 0);

dv1 = v1t' - v1i;
dv2 = v2f - v2t';
dvtot = norm(dv1) + norm(dv2);

end
