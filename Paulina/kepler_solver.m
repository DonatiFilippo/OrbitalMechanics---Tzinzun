function E = kepler_solver(varargin)
    % Inputs:
    % t   - Tiempo actual [s]
    % e   - Excentricidad orbital (0 <= e < 1)
    % a   - Semieje mayor [km]
    % mu  - Parámetro gravitacional [km^3/s^2]
    % t0  - Tiempo inicial [s]
    % E0  - Anomalía excéntrica inicial [rad]
    %
    % Output:
    % E   - Anomalía excéntrica [rad]

    default_tol = 1e-6;

    if nargin == 1
        input_vector = varargin{1}; 
        t = input_vector(1);       
        e = input_vector(2);
        a = input_vector(3);   
        mu = input_vector(4);
        t0 = input_vector(5);
        E0 = input_vector(6);

        if length(input_vector) == 7
            root_method = input_vector(7);
        else
            root_method = 0; % Default value used if input not given
        end

        % Set default tolerance if parameter not provided
         if length(input_vector) == 8
            tol = input_vector(8);
        else
            tol = default_tol; % Default value used if input not given
         end


    elseif nargin >= 6
        t = varargin{1};       
        e = varargin{2};
        a = varargin{3};   
        mu = varargin{4};
        t0 = varargin{5};
        E0 = varargin{6};

        if nargin == 7
            root_method = varargin{7};
        else
            root_method = 0;
        end

         if nargin == 8
            tol = varargin{8};
        else
            tol = default_tol; % Valor predeterminado
        end
    end

    % Define the period
    T = 2*pi*sqrt( a^3/mu ); % Orbital period [s] - For n_ptb orbit 

    % Define to obtain elapsed time
    delta_t = t - t0;

    % Define full orbits
    k = floor(delta_t / T);

    % Define t_bar (remaining time)
    t_bar = delta_t - (k*T);

    % Calcular el movimiento medio (n)
    n = sqrt(mu / a^3); % Velocidad media angular [rad/s]

    % Calcular la anomalía media (M)
    M = n * t_bar; % [rad]

    % Función implícita para resolver la ecuación de Kepler
    F = @(E) E - e * sin(E) - M;

    % Usar fsolve para resolver
    options = optimset('Display', 'off'); % Opciones para suprimir salida
    EK = fsolve(F, E0, options); % Resolver para E usando una aproximación inicial

    if root_method == 0
        E_bar = EK;


    elseif root_method == 1
        % Using fzero to solve the Eq.
        options_fz = optimset('TolX', tol);
        E_bar = fzero(EK,E0,options_fz);
        % fprintf('Solution using fzero:  %.8f\n', E_bar); % Output results 

    elseif root_method == 2
        % Using fsolve to solve the Eq.
        options_fs = optimset('Display', 'none','TolFun', tol);
        E_bar = fsolve(EK, E0, options_fs);
        % fprintf('Solution using fsolve: %.8f\n', E_bar);% Output results
    end
    
    % Add full orbits contributions
    E_tot = k*2*pi + E_bar; % Total eccentric anomaly

    % Return the result
    E = E_tot;  % Return fzero result as default



end

