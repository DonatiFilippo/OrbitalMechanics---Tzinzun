function [s_dot] = twoBodyProblemPert(~, s, mu, J, R)
%TWOBODYPROBLEMPERT Summary of this function goes here
%   Detailed explanation goes here
rr = s(1:3);
vv = s(4:6);

r = norm(rr);
aa = (3/2) * ((J*mu*R^2)/r^4) * [((rr(1)/r) * ((5*rr(3)^2 / r^2) - 1));
                                ((rr(2)/r) * ((5*rr(3)^2 / r^2) - 1));
                                ((rr(3)/r) * ((5*rr(3)^2 / r^2) - 3))];

s_dot = vertcat(vv, ((-mu/r^3)*rr) + aa);
end

