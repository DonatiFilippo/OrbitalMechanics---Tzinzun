function dkep = Gauss_TNH(t, kep, mu, R, J2, muM, date)
%UEQUATION OF MOTION FOR KEPLERIAN PARAMETERS PROPAGATION
% kep = [a, e, i, OM, om, theta]
% a perturbing acceleration in TNH fram

a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
theta = kep(6);

% Evaluate perturbing acceleration
[rp, vp] =  parorb2rv (kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), mu);
ap_ECI = aJ2(rp, mu, R, J2) + a3B_M(t, rp, muM, date);
ap = GE2TNH(rp, vp) * ap_ECI;

% Evaluate equations of motion (Keplerian parameters evolution)
p = a*(1-e^2);
r = p / ( 1+e*cos(theta));
v = sqrt(2*mu/r - mu/a);
n = sqrt(mu/(a^3));
b = a*sqrt(1-e^2);
h = n * a * b;

dkep = [2*(a^2)*v*ap(1,1)/mu;
        1/v * ( 2 * (e + cos(theta)) * ap(1,1) - r/a * sin(theta) * ap(2,1) );
        r * cos(theta + om)/h * ap(3,1);
        r * sin(theta + om)/(h * sin(i)) * ap(3,1);
        1/(e*v) * ( 2*sin(theta)*ap(1,1) + (2*e + r/a*cos(theta))*ap(2,1) ) - (r*sin(theta+om)*cos(i))/(h*sin(i))*ap(3,1);
        h/r^2 - 1/(e*v) * (2 * sin(theta) * ap(1,1) + (2*e + r/a * cos(theta))*ap(2,1)) ];
end