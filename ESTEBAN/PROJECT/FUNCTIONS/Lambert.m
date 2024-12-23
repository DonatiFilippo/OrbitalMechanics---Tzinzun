function dvtot = vtot (t1, t2)

% Initial condition

[kep1,~] = uplanet(t1, 3);
[kep2,~] = uplanet(t2, 4);

mu = astroConstants(4);

[r1, v1i] = kep2cart(kep1(1), kep1(2), kep1(3), kep1(4), kep1(5), kep1(6), mu);
[r2, v2f] = kep2cart(kep2(1), kep2(2), kep2(3), kep2(4), kep2(5), kep2(6),mu);

% Set time quantities
dt = (t2 - t1)*24*60*60;

% Solver
[~, ~, ~, ~, v1t, v2t, ~, ~] = lambertMR(r1, r2, dt, mu, 0, 0, 0, 0);

dv1 = v1t' - v1i;
dv2 = v2f - v2t';
dvtot = norm(dv1) + norm(dv2);

end