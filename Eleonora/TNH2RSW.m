function A_RSW = TNH2RSW(kep, mu)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

a = kep(1);
e = kep(2);
theta = kep(6);

p = a*(1-e^2);
r = p / ( 1+e*cos(theta));
v = sqrt(2*mu/r - mu/a);
n = sqrt(mu/(a^3));
b = a*sqrt(1-e^2);
h = n * a * b;

A_RSW = h/(p*v)*[e*sin(theta), -(1+e*cos(theta)), 0;
                  1+e*cos(theta), e*sin(theta), 0;
                  0, 0, (p*v)/h];
end