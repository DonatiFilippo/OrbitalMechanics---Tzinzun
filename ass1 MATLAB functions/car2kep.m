function [a, e, i, Omega, omega, nu] = car2kep(varargin)
%% FUNCTION DESCRIPTION
%
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function converts a state vector (position and velocity in Cartesian 
% coordinates) into Keplerian orbital elements. It is applicable for any orbit 
% around a central body, assuming a Newtonian two-body problem.
%
%--------------------------------------------------------------------------
% INPUTS:
%   varargin    [1x6] or [1x2]
%       - A single input vector containing position and velocity components:
%         [r_x, r_y, r_z, v_x, v_y, v_z] [km, km/s].
%       - Alternatively, the position and velocity vectors can be provided 
%         separately as [r] and [v] ([km] and [km/s]).
%
%--------------------------------------------------------------------------
% OUTPUTS:
%   a       [1x1]   Semi-major axis                                [km]
%   e       [1x1]   Eccentricity                                   [-]
%   i       [1x1]   Inclination                                    [deg]
%   Omega   [1x1]   Right ascension of the ascending node (RAAN)   [deg]
%   omega   [1x1]   Argument of periapsis                          [deg]
%   nu      [1x1]   True anomaly                                   [deg]
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

    if nargin == 1 && isvector(varargin{1})
        State = varargin{1};
        r = State(1:3);
        v = State(4:6);
    elseif nargin == 2
        r = varargin{1};
        v = varargin{2};
    end
    
    % 1. Normas de los vectores de posición y velocidad
    r_norm = norm(r);
    v_norm = norm(v);
    
    % 2. Vector de momento angular especifico (h)
    h = cross(r, v);
    h_norm = norm(h);
    
    % 3. Cálculo de la inclinación (i)
    i = acos(h(3) / h_norm);
    
    % 4. Nodo ascendente
    N = cross([0, 0, 1], h);
    N_norm = norm(N);
    
    if N_norm ~= 0
        Omega = acos(N(1) / N_norm);
        if N(2) < 0
            Omega = 2 * pi - Omega;
        end
    else
        Omega = 0;
    end
    
    % 5. Vector de excentricidad
    e_vec = (1 / 398600) * ((v_norm^2 - 398600 / r_norm) * r - dot(r, v) * v);
    e = norm(e_vec);
    
    % 6. Argumento del periapsis (omega)
    if N_norm ~= 0
        omega = acos(dot(N, e_vec) / (N_norm * e));
        if e_vec(3) < 0
            omega = 2 * pi - omega;
        end
    else
        omega = 0;
    end
    
    % 7. Análisis verdadero (nu)
    nu = acos(dot(e_vec, r) / (e * r_norm));
    if dot(r, v) < 0
        nu = 2 * pi - nu;
    end
    
    % 8. Semi-eje mayor (a)
    epsilon = (v_norm^2) / 2 - 398600 / r_norm; % Energía específica
    if e ~= 1 % Para órbitas elípticas
        a = -398600 / (2 * epsilon);
    else % Para órbitas parabólicas
        a = Inf;
    end
    
    % Convertir los ángulos de radianes a grados
    i = rad2deg(i);
    Omega = rad2deg(Omega);
    omega = rad2deg(omega);
    nu = rad2deg(nu);
end
