function vrot = Rotate(v, u ,delta)

vrot = v*cos(delta) + cross(u,v)*sin(delta) + dot(u,v)*u*(1-cos(delta));

end