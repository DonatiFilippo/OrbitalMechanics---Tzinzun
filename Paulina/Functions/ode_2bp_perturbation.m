function dy = ode_2bp_perturbation(varargin)
% ode_2bp ODE system for the two-body problem (Keplerian motion)
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

if nargin == 1
        input_vector = varargin{1}; 
        t = input_vector(1);       
        y = input_vector(2);
        mu = input_vector(3);   
        J2 = input_vector(4);
        Re = input_vector(5);
    elseif nargin == 5
        t = varargin{1};
        y = varargin{2}; 
        mu = varargin{3}; 
        J2 = varargin{4};
        Re = varargin{5};
    end


% Position and velocity
r = y(1:3);
v = y(4:6);
% Distance from the primary
rnorm = norm(r);

% Perturbation
x = r(1); 
y2 = r(2); 
z = r(3);

fact = (3/2) * J2*(mu*Re^2)/rnorm^4;
aj_i =  fact * (x / rnorm) * (5*(z^2/rnorm^2) - 1);
aj_j =  fact * (y2 / rnorm) * (5*(z^2/rnorm^2) - 1);
aj_k =  fact * (z / rnorm) * (5*(z^2/rnorm^2) - 3);
aj = [aj_i; aj_j; aj_k];

% Set the derivatives of the state
dy = [ v; (-mu/rnorm^3)*r + aj] ;
end
