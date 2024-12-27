function dvtot = Lambert (t1, t2, planet1, planet2)
 muS = astroConstants(4);

% Initial condition

if planet1 <= 10
    [kep1,~] = uplanet(t1, planet1);
else
    [kep1,~] = ephAsteroids(t1, planet1);
end

if planet2 <= 10
    [kep2,~] = uplanet(t2, planet2);
else
    [kep2,~] = ephAsteroids(t2, planet2);
end

muS = astroConstants(4);

[r1, v1i] = kep2cart(kep1, muS);
[r2, v2f] = kep2cart(kep2,muS);

% Set time quantities
dt = (t2 - t1)*24*60*60;

% Solver
[~, ~, ~, ~, v1t, v2t, ~, ~] = lambertMR(r1, r2, dt, muS, 0, 0, 0, 0);

dv1 = v1t' - v1i;
dv2 = v2f - v2t';
dvtot = norm(dv1) + norm(dv2);

end