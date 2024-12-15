function [a] = RepeatGroundTrack(k, m, mu, omega_E)
% computeRepeatingGroundTrack.m - Computes the semimajor axis for a repeating ground track.
%
% Inputs:
%   k         - Number of satellite revolutions
%   m         - Number of Earth rotations
%   mu        - Gravitational parameter of Earth (m^3/s^2)
%   omega_E   - Earth's angular velocity (rad/s)
%   t0        - Initial time (s) (not explicitly used here, kept for structure consistency)
%
% Output:
%   a         - Semimajor axis of the orbit (m)
%
% Author:
%   Maria Paulina Pantoja Gavidia
% ------------------------------------------------------------------------

% Step 1: Earth's sidereal rotation period [s]
T_E = 2 * pi / omega_E;

% Step 2: Orbital period for a repeating ground track
T = T_E * (m / k);

% Step 3: Compute semimajor axis using Kepler's Third Law
a = (mu * T^2 / (4 * pi^2))^(1/3);

end
