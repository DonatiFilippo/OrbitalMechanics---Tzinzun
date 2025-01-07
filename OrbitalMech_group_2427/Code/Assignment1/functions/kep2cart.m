function [r_vect,v_vect] = kep2cart(varargin)
%
% Converter: Keplerian Parameters ---> Cartesian Coordinates & Velocities
%
%DESCRIPTION:
%This code provides a conversion from the Keplerian Parameters of the orbit
%to cartesian coordinates vector (r_vect) and velocities (v_vect).
%
%--------------------------------------------------------------------------
% INPUTS:
%   a          [1x1]       Semi-Major Axis           [km]
%   e          [1x1]       Eccentricity              [-]
%   i          [1x1]       Inclination               [rad]
%   OM         [1x1]       RAAN                      [rad]
%   om         [1x1]       Argument of Periapsis     [rad]
%   th         [1x1]       True Anomaly              [rad]
%   mu         [1x1]       Planetary Constant        [km^3][sec-2]
%--------------------------------------------------------------------------
% OUTPUTS:
%   r_vect     [3x1]       Position Vector           [km]
%   v_vect     [3x1]       Velocity Vector           [km/s]
%--------------------------------------------------------------------------
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
if nargin == 1
    % Formato: kep2cart(kep)
    if isvector(varargin{1})
        Kep = varargin{1};
        a = Kep(1);
        e = Kep(2);
        i = Kep(3);
        OM = Kep(4);
        om = Kep(5);
        th = Kep(6);
        mu = 398600; % Valor por defecto para la Tierra [km^3/s^2]
    end

elseif nargin == 2
    % Formato: kep2cart(kep, mu)
    if isvector(varargin{1}) && isscalar(varargin{2})
        Kep = varargin{1};
        a = Kep(1);
        e = Kep(2);
        i = Kep(3);
        OM = Kep(4);
        om = Kep(5);
        th = Kep(6);
        mu = varargin{2};
    end

elseif nargin == 7
    % Formato: kep2cart(kep1(1), kep1(2), ..., kep1(6), mu)
    a = varargin{1};
    e = varargin{2};
    i = varargin{3};
    OM = varargin{4};
    om = varargin{5};
    th = varargin{6};
    mu = varargin{7};
    end


%% Conversion Routine
%Semi-Latus Rectum
p = a*(1-e.^2);
%Radial Distance
r=p/(1+e*cos(th));
%Position and Velocity Vectors in Perifocal Coordinates
r_pf = r*[cos(th);sin(th);0];
v_pf = sqrt(mu/p)*[-sin(th);e+cos(th);0];
%Rotation Matrix
R = [cos(om)*cos(OM)-sin(om)*cos(i)*sin(OM) -sin(om)*cos(OM)-cos(om)*cos(i)*sin(OM) sin(i)*sin(OM);
     cos(om)*sin(OM)+sin(om)*cos(i)*cos(OM) -sin(om)*sin(OM)+cos(om)*cos(i)*cos(OM) -sin(i)*cos(OM);
     sin(om)*sin(i)                         cos(om)*sin(i)                          cos(i)];
%Position and Velocity Vectors in Cartesian Coordinates
r_vect = R*r_pf;
v_vect = R*v_pf;
end
