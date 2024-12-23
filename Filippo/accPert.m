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
r = p / ( 1+e*cos(f));
h = sqrt(p*mu);
v = sqrt(2*mu/r - mu/a);

aa = (3/2) * ((J*mu*R^2)/r^4) * [((rr(1)/r) * ((5*rr(3)^2 / r^2) - 1));
                                ((rr(2)/r) * ((5*rr(3)^2 / r^2) - 1));
                                ((rr(3)/r) * ((5*rr(3)^2 / r^2) - 3))];
amcar = accMoon(t, s, muM, mu,start) + aa;
% amtnh = car2TNH(rr, vv)* amcar;
% amrs = [(h/(p*v))*[e*sin(f), -(1 + e*sin(f)), 0]; 
%         (h/(p*v))*[1+e*cos(f), e*sin(f), 0];
%         [0, 0, 1]];
% am = amrs * amtnh;

R3_thetaw = [cos(om+f) sin(om+f) 0; -sin(om+f) cos(om+f) 0; 0 0 1];
R1_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R3_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
am = R3_thetaw*R1_i*R3_OM*amcar;
accp = am;
end

