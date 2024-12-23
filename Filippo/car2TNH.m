function A_TNH = car2TNH(r,v)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
t_uv = v / norm(v);
h = cross(r, v);
h_uv = h / norm(h);
n_uv = cross(h_uv, t_uv);

A_TNH = [t_uv, n_uv, h_uv]';
end