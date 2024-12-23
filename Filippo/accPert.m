function [accp] = accPert(t, s, mu, J, R, muM,start)
%ACCPERT Summary of this function goes here
%   Detailed explanation goes here
a = s(1);
e = s(2);
i = s(3);
OM = s(4);
om = s(5);
f = s(6);

[rr, vv] = kep2car(a,e, i, OM, om, f, mu);
p = a*(1-(e^2));
r = p/(1+e*cos(f));
h = sqrt(p*mu);
v = sqrt((2*mu/r) - mu/a);

aa = -(3/2) * ((J*mu*(R^2))/(r^4)) * [(1 - 3 * ((sin(i))^2)) * ((sin(f+om))^2);
                                           (((sin(i))^2) * (sin(2*(f+om))));
                                           (sin(2*i) * sin(f+om))];
amcar = accMoon(t, s, muM, mu,start);
amtnh = car2TNH(rr, vv) * amcar;
amrs = [h/(p*v)*[e*sin(f), -(1 + e*sin(f)), 0]; 
        h/(p*v)*[1+e*cos(f), e*sin(f), 0];
        [0, 0, 1]];
am = amrs * amtnh;
accp = aa + am;
end

