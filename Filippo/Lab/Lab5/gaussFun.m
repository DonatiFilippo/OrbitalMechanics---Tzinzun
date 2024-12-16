function [s_dot] = gaussFun(t, s, acc, mu)
%GAUSSFUN Summary of this function goes here
%   Detailed explanation goes here

a = s(1);
e = s(2);
i = s(3);
OM = s(4);
w = s(5);
theta = s(6);
th = theta;
om = w;

p = a*(1-e^2);
r = p/(1+e*cos(theta));
h = sqrt(p*mu);

ar = acc(1);
as = acc(2);
aw = acc(3);

a_dot = (2*a^2/h)*(e*sin(th)*ar+p*as/r);
e_dot = (1/h)*(p*sin(th)*ar+((p+r)*cos(th)+r*e)*as);
i_dot = r*cos(th+om)*aw/h;
OM_dot = r*sin(th+om)*aw/(h*sin(i));
om_dot = (1/(h*e))*(-p*cos(th)*ar+(p+r)*sin(th)*as)-r*sin(th+om)*cos(i)*aw/(h*sin(i));
th_dot = h/r^2+(1/(e*h))*(p*cos(th)*ar-(p+r)*sin(th)*as);

s_dot = vertcat(a_dot, e_dot, i_dot, OM_dot, om_dot, th_dot);
end

