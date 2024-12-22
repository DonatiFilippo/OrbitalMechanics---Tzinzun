function [acc_vec] = acc_pert(~, s, J, mu, R)
%ACC_PERT Summary of this function goes here
%   Detailed explanation goes here

a = s(1);
e = s(2);
i = s(3);
om = s(5);
f = s(6);

p = a*(1-(e^2));
r = p/(1+e*cos(f));

acc_vec = -(3/2) * ((J*mu*(R^2))/(r^4)) * [(1 - 3 * ((sin(i))^2)) * ((sin(f+om))^2);
                                           (((sin(i))^2) * (sin(2*(f+om))));
                                           (sin(2*i) * sin(f+om))];
                                    
end

