function dy = ode_2bp(varargin)
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
%
% -------------------------------------------------------------------------

% Check if entrance is vector o dependly values
 if nargin == 1
        input_vector = varargin{1}; 
        t = input_vector(1);       
        y = input_vector(2:end-1);
        mu = input_vector(end);   
    elseif nargin == 3
        t = varargin{1};
        y = varargin{2}; 
        mu = varargin{3}; 
    end

% Position and velocity
r = y(1:3);
v = y(4:6);

% Distance from the primary
rnorm = norm(r);
% Set the derivatives of the state
dy = [ v
 (-mu/rnorm^3)* r ];
end


