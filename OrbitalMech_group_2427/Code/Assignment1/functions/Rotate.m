function vrot = Rotate(v, u ,delta)
%% FUNCTION DESCRIPTION
%
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function performs a rotation of a vector `v` around a unit vector `u` 
% by an angle `delta` (in radians). It uses the Rodrigues' rotation formula 
% for efficient computation of the rotated vector.
%
%--------------------------------------------------------------------------
% INPUTS:
%   v          [3x1]       Input vector to be rotated                 [km/s]
%   u          [3x1]       Unit vector defining the axis of rotation  [-]
%   delta      [1x1]       Rotation angle (in radians)                [rad]
%
%--------------------------------------------------------------------------
% OUTPUTS:
%   vrot       [3x1]       Rotated vector                             [km/s]
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


vrot = v*cos(delta) + cross(u,v)*sin(delta) + dot(u,v)*u*(1-cos(delta));

end
