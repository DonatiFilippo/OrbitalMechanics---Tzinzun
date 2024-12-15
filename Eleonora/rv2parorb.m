function [a, e, i, OM, om, theta] = rv2parorb(rr, vv, mu)

%Transformation from Cartesian state to orbital elements
%
%Input arguments:
%rr             [3x1]  position vector              [km]
%vv             [3x1]  velocity vector              [km/s]
%mu             [1x1]  gravitational parameter      [km^3/s^2]
%
%-------------------------------------------------------------------------
%Output arguments:
%a              [1x1]  semi-major axis              [km]
%e              [1x1]  eccentricity                 [-]
%i              [1x1]  inclination                  [rad]
%OM             [1x1]  RAAN                         [rad]
%om             [1x1]  pericenter anomaly           [rad]
%theta          [1x1]  true anomaly                 [rad]
%-------------------------------------------------------------------------

%versori del sistema geocentrico
I=[1; 0; 0];
J=[0; 1; 0];
K=[0; 0; 1];

%moduli dei vettori posizione (rr) e velocità (vv)
r=norm(rr);
v=norm(vv);

v_r=dot(rr,vv)/r; %velocità radiale

%emiasse maggiore a dalla formula dell'energia
a=1/(2/r-(v^2)/mu);

hh=cross(rr,vv); %vettore momento angolare
h=norm(hh); %modulo del vettore momento angolare
i=acos(dot(hh,K)/h); %inclinazione

%se l'orbita è equatoriale poniamo NN=I, l'ascensione retta del nodo
%ascendente risulterà nulla
if (i == 0) 
    NN=I;
else
    NN=cross(K,hh); %linea dei nodi
end
N=norm(NN); %magnitude of node vector

%ascensione retta del nodo ascendente
if(dot(NN,J)>=0)
    OM=acos(dot(NN,I)/N);
else
    OM=2*pi-acos(dot(NN,I)/N);
end

%eccentricità
ee=(cross(vv,hh))/mu-(rr/r);

%se l'orbita è circolare (e=0) poniamo il vettore eccentricità coincidente
%con la linea dei nodi
if (norm(ee) <= 1e-8)
    ee=NN;
end
e=norm(ee);

%argomento del pericentro
if(dot(ee,K)>=0)
    om=acos(dot(NN,ee)/(N*e));
else
    om=2*pi-acos(dot(NN,ee)/(N*e));
end

%anomalia vera
if(v_r>=0)
    theta=acos(dot(ee,rr)/(e*r));
else
    theta=2*pi-acos(dot(ee,rr)/(e*r));
end
   
end