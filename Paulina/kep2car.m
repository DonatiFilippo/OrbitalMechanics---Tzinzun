function [r0, v0] = kep2car(varargin)

    % kep2car - Converts Keplerian elements to Cartesian state vectors.
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
        Kep = varargin{1};
        a = Kep(1);
        e = Kep(2);
        i = Kep(3);
        Omega = Kep(4);
        omega = Kep(5);
        nu = Kep(6);
        mu = Kep(7);
    elseif nargin == 7
        a = varargin{1};
        e = varargin{2};
        i = varargin{3};
        Omega = varargin{4};
        omega = varargin{5};
        nu = varargin{6};
        mu = varargin{7};
    end

    % Convert angles from degrees to radians
    i = deg2rad(i);
    Omega = deg2rad(Omega);
    omega = deg2rad(omega);
    nu = deg2rad(nu);
    
    % Calculate the semi-latus rectum
    p = a * (1 - e^2);
    
    % Position magnitude in the perifocal coordinate system
    r = p / (1 + e * cos(nu));
    
    % Position and velocity in perifocal coordinates
    r_perifocal = [r * cos(nu); r * sin(nu); 0];
    v_perifocal = sqrt(mu / p) * [-sin(nu); e + cos(nu); 0];
    
    % Construct each element of the rotation matrix cR
    R3_W = [ cos(Omega),  sin(Omega),  0;
            -sin(Omega),  cos(Omega),  0;
             0,          0,         1];
         
    R1_i = [1,  0,          0;
            0,  cos(i),   sin(i);
            0, -sin(i),   cos(i)];
         
    R3_w = [ cos(omega),  sin(omega),  0;
            -sin(omega),  cos(omega),  0;
             0,       0,       1];
         
    % Matriz de rotación completa del sistema PQW al sistema ECI
    Q_pqw2eci = R3_W * R1_i * R3_w;
    
    % 3. Transformar la posición y la velocidad de PQW a ECI
    r0 = Q_pqw2eci * r_perifocal; % Vector de posición en el sistema ECI (km)
    v0 = Q_pqw2eci * v_perifocal; % Vector de velocidad en el sistema ECI (km/s)
end