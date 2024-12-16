function dkep = eq_motion( t,kep, acc_pert_fun, J2, mu, Re )
% Evaluate the perturbing accelerations

%s = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6));
% acc_pert_vec = acc_pert_fun( t,kep, J2, mu, Re );
acc_pert_vec = acc_pert_fun;


% Evaluate the equations of motion (Cartesian or Keplerian),
% function of t, s, acc_pert_vec, and parameters

% Acceleration from inertial frame to RSW frame
a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
w = kep(5);
theta = kep(6);
ar = acc_pert_vec(1);
as = acc_pert_vec(2);
aw = acc_pert_vec(3);


% Gauss planetary equations RSW frame
r = a*(1-e^2)/(1+e*cos(theta));
% v = sqrt(2.*mu/r - mu/a);
b = a*sqrt(1-e^2);
p = a*(1-e^2);
n = sqrt(mu/a^3);
h = n*a*b;

da = (2*a^2/h)*(e*sin(theta)*ar + p*as/r);
de = (1/h)*(p*sin(theta)*ar + ((p+r)*cos(theta)+r*e)*as);
di = r*cos(theta+w)*aw/h;
dOM = r*sin(theta+w)*aw/(h*sin(i));
dw = (1/(h*e))*(-p*cos(theta)*ar + (p+r)*sin(theta)*as)-...
    r*sin(theta+w)*cos(i)*aw/(h*sin(i));
dtheta = h/r^2 + (1/(e*h))*(p*cos(theta)*ar-(p+r)*sin(theta)*as);

dkep = [da de di dOM dw dtheta]'; 

display(t)

end