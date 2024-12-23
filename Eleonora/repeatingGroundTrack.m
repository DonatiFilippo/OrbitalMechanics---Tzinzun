function [a] = repeatingGroundTrack(m, k, wE, mu)
% repeatingGroundTrack computes the required semi-major axis a for a
% repeating ground track given m and k
%
% PROTOTYPE: a = repeatingGroundTrack(m, k, wE, mu)
%
% INPUT: 
% m[1x1]   Primary's revolutions                              [-]
% k[1x1]   Satellite's revolutions                            [-]
% wE[1x1]  Primary's angular velocity                         [rad/T]
% mu[1x1]  Gravitational parameter of the primary             [L^3/T2]
%
% OUTPUT:
% a[1x1]   Semi-major axis to have repeating GT                [L]
%
% CONTRIBUTORS:
%  Eleonora Domenichelli
%
% VERSION:
%  2024/11/17: First version
%
%-------------------------------------------------------------------------
a = (mu/(wE*(k/m))^2)^(1/3);
end