function dkep = Gauss_RSW(t, kep, mu, R, J2, muM, date)
%UEQUATION OF MOTION FOR KEPLERIAN PARAMETERS PROPAGATION
% kep = [a, e, i, OM, om, theta]
% a perturbing acceleration in TNH fram

a = kep(1);
e = kep(2);
i = kep(3);
om = kep(5);
theta = kep(6);

% Evaluate perturbing acceleration
[rp, vp] =  parorb2rv (kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), mu);
ap_ECI = aJ2(rp, mu, R, J2) + a3B_M(t, rp, muM, date);
ap_TNH = car2TNH(rp, vp) * ap_ECI;
ap = TNH2RSW(kep, mu) * ap_TNH;

% Evaluate equations of motion (Keplerian parameters evolution)
p = a*(1-e^2);
r = p / ( 1+e*cos(theta));
n = sqrt(mu/(a^3));
b = a*sqrt(1-e^2);
h = n * a * b;

dkep = [2*(a^2)/h * (e*sin(theta)*ap(1,1) + p/r * ap(2,1));
        1/h * (p*sin(theta)*ap(1,1) + ((p+r)*cos(theta) + r*e)*ap(2,1));
        r*cos(theta+om)/h * ap(3,1);
        r*sin(theta+om)/(h*sin(i)) * ap(3,1);
        1/(h*e) * (-p*cos(theta)*ap(1,1) + (p+r)*sin(theta)*ap(2,1)) - (r*sin(theta+om)*cos(i))/(h*sin(i)) * ap(3,1);
        h/(r^2) + 1/(e*h) * (p*cos(theta)*ap(1,1) - (p+r)*sin(theta)*ap(2,1))];
end