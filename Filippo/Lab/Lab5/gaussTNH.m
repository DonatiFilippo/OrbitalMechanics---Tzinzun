function [s_dot] = gaussTNH(t, s, acc, mu)
%GAUSSTNH Summary of this function goes here
%   Detailed explanation goes here

a = s(1);
e = s(2);
i = s(3);
om = s(5);
th = s(6);

p = a*(1-(e^2));
r = p / ( 1+e*cos(th));
h = sqrt(p*mu);
v = sqrt(2*mu/r - mu/a);



at = acc(1);
an = acc(2);
ah = acc(3);

a_dot = (2*(a^2)*v*at/mu);
e_dot = 1/v * (2*(e+cos(th))*at - r*sin(th)*an/a);
i_dot = r*cos(th+om)*ah/h;
OM_dot = r*sin(th+om)*ah/(h*sin(i));
om_dot = (1/(e*v))*(2*sin(th)*at + (2*e+r*cos(th)/a)*an) - r*sin(th+om)*cos(i)*ah/(h*sin(i));
th_dot = h/(r^2)-(1/(e*v))*(2*sin(th)*at+(2*e+r*cos(th)/a)*an);



s_dot = vertcat(a_dot, e_dot, i_dot, OM_dot, om_dot, th_dot);

end

