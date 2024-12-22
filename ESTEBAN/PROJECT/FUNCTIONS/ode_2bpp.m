function dy = ode_2bpp( ~ , y , mu , J2 , Re )

r = y(1 : 3);
v = y (4 : 6);

rnorm = norm(r);

aj2 = 1.5*(J2*mu*Re^2/rnorm^4)*[(y(1)/rnorm)*(5*(y(3)^2/rnorm^2)-1);
    (y(2)/rnorm)*(5*(y(3)^2/rnorm^2)-1);
    (y(3)/rnorm)*(5*(y(3)^2/rnorm^2)-3)];

dy = [v 
        (-mu/rnorm^3)*r + aj2];

end