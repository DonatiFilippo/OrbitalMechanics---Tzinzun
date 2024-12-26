function dvtot = Lambert (t1, t2, planet1, planet2)

% Initial condition

if planet1 <= 10
    [kep1,~] = uplanet(t1, planet1);
else
    [kep1,~] = ephAsteroids(t1, planet1);
end

if planet2 <= 10
    [kep2,~] = uplanet(t1, planet1);
else
    [kep2,~] = ephAsteroids(t1, planet1);
end

muS = astroConstants(4);

[r1, v1i] = kep2cart(kep1(1), kep1(2), kep1(3), kep1(4), kep1(5), kep1(6), muS);
[r2, v2f] = kep2cart(kep2(1), kep2(2), kep2(3), kep2(4), kep2(5), kep2(6),muS);

% Set time quantities
dt = (t2 - t1)*24*60*60;

% Solver
[~, ~, ~, ~, v1t, v2t, ~, ~] = lambertMR(r1, r2, dt, muS, 0, 0, 0, 0);

dv1 = v1t' - v1i;
dv2 = v2f - v2t';
dvtot = norm(dv1) + norm(dv2);

end