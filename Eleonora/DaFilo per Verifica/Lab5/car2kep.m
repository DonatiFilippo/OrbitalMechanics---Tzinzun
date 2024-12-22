function [a, e, i, OM, w, f] = car2kep(rr, vv, mu)
%CAR2KEP Summary of this function goes here
%   Detailed explanation goes here
r = norm(rr);
v = norm(vv);
vr = dot(vv*rr)/r;
hh = cross(r,v);
h = norm(hh);
zz = [0,0,1];

i = acos(hh(3)/h);
NN = cross(zz, hh/h);

if(NN(2)>= 0)
    OM = acos(NN(1));
else
    OM = 2*pi - acos(NN(1));
end

if(i == 0)
    OM = 0;
    NN = [1, 0, 0];
end

ee = (1/mu) * (((v^2)-mu/r)*rr -vr*r*vv);
e = norm(ee);

if(ee == 0)
    w = 0;
elseif(ee(3) >= 0)
    w = acos(dot(NN,ee)/norm(NN)*e);
else
    w = 2*pi - acos(dot(NN,ee)/norm(NN)*e);
end

if(vr >0)
    f = acos(dot(ee,rr)/(e*r));
elseif(vr <0)
    f = 2*pi - acos(dot(ee,rr)/(e*r));
else
    if(rr == ee)
        f = 0;
    else
        f = pi;
    end
end

a = (-mu/2)/(v^2/2 - mu/r);

end

