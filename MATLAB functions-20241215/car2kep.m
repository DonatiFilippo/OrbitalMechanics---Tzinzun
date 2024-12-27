function [a, e, i, Omega, omega, nu] = car2kep(varargin)
    % CAR2KEP - Convierte vectores de estado cartesianos a elementos keplerianos.
    %
    % Inputs:
    %   a       [1]   Semi-major axis                      [km]
    %   e       [1]   Eccentricity                         [-]
    %   i       [1]   Inclination                          [deg]
    %   Omega   [1]   Right ascension of the ascending node [deg]
    %   omega   [1]   Argument of perigee                  [deg]
    %   nu       [1]   True anomaly                         [deg]
    %   mu      [1]   Gravitational parameter of the primary body [km^3/s^2]
    % Input can be a vector
    %
    % OUTPUT:
    %   r       [3x1]   Position vector in Cartesian coordinates [km]
    %   v       [3x1]   Velocity vector in Cartesian coordinates [km/s]
    %
    %
    % Author
    %  Maria Paulina Pantoja Gavidia
    % -------------------------------------------------------------------------
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
