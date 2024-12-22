function dy = ode_2bp(~, y, mu)
%ODE_2BP Summary of this function goes here
%   Detailed explanation goes here
rr = y(1:3);
vv = y(4:6);

r = norm(rr);

dy = [vv;
        (-mu/r^3)*rr];
end

