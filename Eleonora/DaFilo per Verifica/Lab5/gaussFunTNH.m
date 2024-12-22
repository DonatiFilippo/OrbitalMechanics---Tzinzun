function [s_dot] = gaussFunTNH(~, s, acc, mu)
%GAUSSFUNTNH Summary of this function goes here
%   Detailed explanation goes here outputArg1 = inputArg1;

a = s(1);
e = s(2);
i = s(3);
w = s(5);
f = s(6);


p = a*(1-e^2);
r = p/(1+e*cos(theta));
h = sqrt(p*mu);
v = sqrt((2*mu/r)-mu/a);

at = acc(1);
an = acc(2);
ah = acc(3);

a_dot = 2*(a^2)*v*at/mu;
e_dot = 1/v * (2*(e+cos(f))*at - r/a*sin(f)*an);
i_dot = r*cos(f+w)*ah / h;
OM_dot = r*sin(f+w)*ah / (h*sin(i));
w_dot = (1/(e*v))*(2*sin(f)*at + (2*e + r/a*cos(f))*an) - r*sin(f+w)*cos(i)*ah/(h*sin(i));
f_dot = h/(r^2) - (1/(e*v))*(2*sin(f)*at + (2*e+r/a*cos(f))*an);

s_dot = vertcat(a_dot, e_dot, i_dot, OM_dot, w_dot, f_dot);
end

