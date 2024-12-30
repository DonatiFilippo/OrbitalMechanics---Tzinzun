function a = a_repeatingGT(m, k, wE, mu)

% a_repeatingGT - Semi-major axis a for the repeating ground track
% given m and k
%
% PROTOTYPE:
%   a = a_repeatingGTm, k, wE, mu)
%
% DESCRIPTION:
%   The function gives the semi-major axis necessary to obtain a repeated
%   ground track, given the parameters k (revolutions of the satellite between
%   repetitions) and m (rotations of the planet between repetitions).
%
% INPUT: 
%   m[1x1]   Primary's rotations                                [-]
%   k[1x1]   Satellite's revolutions                            [-]
%   wE[1x1]  Primary's angular velocity                         [rad/s]
%   mu[1x1]  Gravitational parameter of the primary             [km^3/s^2]
%
% OUTPUT:
%   a[1x1]   Semi-major axis to have repeating GT               [km]
%
% CONTRIBUTORS:
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo 
%   Domenichelli Eleonora
%
%-------------------------------------------------------------------------

a = (mu/(wE*(k/m))^2)^(1/3);
end