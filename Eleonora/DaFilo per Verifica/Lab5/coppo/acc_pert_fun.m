function acc_pert_vec = acc_pert_fun(kep, J2, mu, Re )
% Evaluate the perturbing accelerations in a given
% reference frame (e.g., TNH, RTH, ECI, etc.)

% State
i = kep(3);
OM = kep(4);
w = kep(5);
theta = kep(6);
s = parOrb2rv(kep(1),kep(2),i,kep(4),w,theta,mu);

% J2 perturbation part

x = s(1);
y = s(2);
z = s(3);
nr = norm([x y z]);
 
aJ2 = (1.5*J2*mu*Re^2/nr^4).*[x/nr.*(5.*(z./nr)^2-1) ...
    y/nr.*(5.*(z./nr)^2-1)  z/nr.*(5.*(z./nr)^2-3)]';

R3_thetaw = [cos(w+theta) sin(w+theta) 0; -sin(w+theta) cos(w+theta) 0; 0 0 1];
R1_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R3_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
acc_pert_vec = R3_thetaw*R1_i*R3_OM*aJ2;
 
end
