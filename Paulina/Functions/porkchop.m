function [delta_v1, delta_v2, r1_arc, r2_arc, r3_arc, v_fb, t_fb] = porkchop(dep_dates_1, arr_dates_1, dep_dates_2, arr_dates_2)

% Inicialización de constantes
muS = astroConstants(4);
muE = astroConstants(13); 

% Inicialización de matrices
delta_v1 = nan(length(dep_dates_1), length(arr_dates_1));
delta_v2 = nan(length(dep_dates_2), length(arr_dates_2));

r1 = zeros(length(dep_dates_1), 3);
v1 = zeros(length(dep_dates_1), 3);
r2 = zeros(length(dep_dates_1), 3);
v2 = zeros(length(dep_dates_1), 3);

r2_2 = zeros(length(dep_dates_2), 3); % Matriz para la segunda etapa
v2_2 = zeros(length(dep_dates_2), 3);
r3 = zeros(length(dep_dates_2), 3);
v3 = zeros(length(dep_dates_2), 3);

% Cálculo de posiciones y velocidades para Mercury y Earth (etapa 1)
for i = 1:length(dep_dates_1)
    [kep_dep, ~] = uplanet(dep_dates_1(i), 1); % Mercury
    [r1(i,:), v1(i,:)] = kep2cart(kep_dep, muS);

    [kep_fb, ~] = uplanet(dep_dates_1(i), 3); % Earth
    [r2(i,:), v2(i,:)] = kep2cart(kep_fb, muS);
end

% Cálculo de posiciones y velocidades para Earth y el asteroide (etapa 2)
for i = 1:length(dep_dates_2)
    [kep_fb2, ~] = uplanet(dep_dates_2(i), 3); % Earth
    [r2_2(i,:), v2_2(i,:)] = kep2cart(kep_fb2, muS);

    [kep_arr, ~, ~] = ephAsteroids(dep_dates_2(i), 30); % Asteroid
    [r3(i,:), v3(i,:)] = kep2cart(kep_arr, muS);
end

% Cálculo de delta-v para la primera etapa (Mercury -> Earth)
for i = 1:length(dep_dates_1)
    for j = 1:length(arr_dates_1)
        tof = (arr_dates_1(j) - dep_dates_1(i)) * 86400; % TOF en segundos
        if tof > 0
            [~,~,~,error,vi,vfb,~,~] = lambertMR(r1(i,:), r2(j,:), tof, muS);
            if error == 0
                delta_v1(i, j) = norm(vi - v1(i,:)) + norm(vfb - v2(j,:));
            else
                delta_v1(i, j) = NaN;
            end
        end
    end
end

% Cálculo de delta-v para la segunda etapa (Earth -> Asteroid)
for i = 1:length(dep_dates_2)
    for j = 1:length(arr_dates_2)
        tof2 = (arr_dates_2(j) - dep_dates_2(i)) * 86400; % TOF en segundos
        if tof2 > 0
            [~,~,~,error2,vfb2,vf,~,~] = lambertMR(r2_2(i,:), r3(j,:), tof2, muS);
            if error2 == 0
                delta_v2(i, j) = norm(vfb2 - v2_2(i,:)) + norm(vf - v3(j,:));
            else
                delta_v2(i, j) = NaN;
            end
        end
    end
end

% Encontrar delta-v mínimo para ambas etapas
[~, idx_min1] = min(delta_v1(:), [], 'omitnan');
[row_min1, col_min1] = ind2sub(size(delta_v1), idx_min1);
r1_arc = r1(row_min1,:);
r2_arc = r2(col_min1,:);

[~, idx_min2] = min(delta_v2(:), [], 'omitnan');
[row_min2, col_min2] = ind2sub(size(delta_v2), idx_min2);
r3_arc = r3(col_min2,:);

% Velocidades en el flyby
v_fb = v2(col_min1,:);
t_fb = dep_dates_1(row_min1);

end
