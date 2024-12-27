function [delta_v1, delta_v2] = porkchop(dep_dates_1, arr_dates_1, dep_dates_2, arr_dates_2)
%% FUNCTION DESCRIPTION
%
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function calculates the delta-v for two consecutive transfer arcs in 
% order to compute the porkchop plot contours, which represent the delta-v 
% cost for various  combinations of departure and arrival dates. 
% These plots help to identify optimal transfer windows.
%
%--------------------------------------------------------------------------
% INPUTS:
%   dep_dates_1    [nx1]   Array of departure dates for the first leg [mjd2000]
%   arr_dates_1    [nx1]   Array of arrival dates for the first leg   [mjd2000]
%   dep_dates_2    [nx1]   Array of departure dates for the second leg [mjd2000]
%   arr_dates_2    [nx1]   Array of arrival dates for the second leg   [mjd2000]
%
%--------------------------------------------------------------------------
% OUTPUTS:
%   delta_v1       [nxm]   Matrix of delta-v values for the first transfer leg [km/s]
%   delta_v2       [pxq]   Matrix of delta-v values for the second transfer leg [km/s]
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

% Inicialización de matrices
delta_v1 = nan(length(dep_dates_1), length(arr_dates_1));
delta_v2 = nan(length(dep_dates_2), length(arr_dates_2));

for i = 1:length(dep_dates_1)
    for j = 1:length(arr_dates_1)
        delta_v1(i, j) = Lambert(dep_dates_1(i), arr_dates_1(j), 1, 3);
    end
end

% Cálculo de delta-v para la segunda etapa (Earth -> Asteroid)
for i = 1:length(dep_dates_2)
    for j = 1:length(arr_dates_2)
        delta_v2(i, j) = Lambert(dep_dates_2(i), arr_dates_2(j), 3, 30);
    end
end

end
