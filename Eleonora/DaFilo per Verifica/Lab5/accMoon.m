function [accM] = accMoon(t, s, muMoon, muEarth, startDate)
%ACCMOON Summary of this function goes here
%   Detailed explanation goes here

a = s(1);
e = s(2);
i = s(3);
OM = s(4);
w = s(5);
f = s(6);


rc3b = ephMoon(startDate+t);
rcs = kep2car(a,e,i,OM,w,f,muEarth);
rs3b = rc3b - rcs;
accM = muMoon* ((rs3b/norm(rs3b)^3) - (rc3b/norm(rc3b)^3));
end

