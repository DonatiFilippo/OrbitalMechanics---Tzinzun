function [delta_v1, delta_v2] = porkchop(dep_dates_1, arr_dates_1, dep_dates_2, arr_dates_2)


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
