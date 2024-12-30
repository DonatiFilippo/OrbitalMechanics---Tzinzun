function plot0rbit(a, e, i, OM, om, th0, thf, dth, mu, mod)

% plot0rbit - Orbit representation
%
% PROTOTYPE:
%   plot0rbit(a, e, i, OM, om, th0, thf, dth, mu, mod)
%
% DESCRIPTION:
%   The function gives as output the orbit representation from Keplerian
%   elements at starting point th0, over the true anomaly interval [th0,
%   thf] with step dth. Note that the central planet is not automatically
%   represented by this function.
%
% INPUT:
%   a [1x1]      Semi-major axis                                 [km]
%   e [1x1]      Eccentricity                                    [-]
%   OM [1x1]     Right ascension of ascending node               [rad]
%   om [1x1]     Pericenter anomaly                              [rad]
%   th0 [1x1]    Starting true anomaly                           [rad]
%   thf [1x1]    Final true anomaly                              [rad]
%   dth [1x1]    True anomaly resolution                         [rad]
%   mu [1x1]     Gravitational parameter of the primary          [km^3/s^2]
%   mod [1x1]    Discrimination parameter for different plots    [-]
%
% OUTPUT:
%   Graphic representation of the orbit from th0 to thf in 3D.
%
% CONTRIBUTORS:
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo 
%   Domenichelli Eleonora
%
%-------------------------------------------------------------------------

% True anomaly vector
t = th0 : dth : thf;

% Position vector in Geocentric equatorial reference frame initialization
l = length(t);
r_GE=zeros(3,l);

for k = 1:l
r_GE(:,k) = parorb2rv (a, e, i, OM, om, t(k), mu);
end

% Option 1: if mod not given to the function, only the orbit is generated
% in a continuous colored line
if nargin == 9
plot3(r_GE(1,:), r_GE(2,:), r_GE(3,:),'LineWidth',2)

else
    % Option 2: if mod == 1 the orbit is dashed in black
    if mod == 1 %dashed orbit in black
        plot3(r_GE(1,:), r_GE(2,:), r_GE(3,:),'k--')

    % Option 3: if mod == 2 the orbit is in continuous colored line and
    % with eccentricity vector and starting point represented 
    elseif mod == 2 
        plot3(r_GE(1,:), r_GE(2,:), r_GE(3,:),'LineWidth',2)

        [r_GE(:,1),vv] = parorb2rv (a, e, i, OM, om,th0, mu);
        hh = cross(r_GE(:,1),vv);
        ee = cross(vv,hh)./mu - r_GE(:,1)./norm(r_GE(:,1));
        ev = ee/norm(ee).*(a*(1-norm(ee)));
        plot3(linspace(0,ev(1),100),linspace(0,ev(2),100),linspace(0,ev(3),100),'y','LineWidth',2)
        plot3(ev(1),ev(2),ev(3),'oy','LineWidth', 6)
    
    % Option 4: if mod == 3 the orbit is dashed in black with eccentricity
    % vector and starting point represented
    elseif mod == 3
        plot3(r_GE(1,:), r_GE(2,:), r_GE(3,:),'k--')

        [r_GE(:,1),vv] = parorb2rv (a, e, i, OM, om,th0, mu);
        hh = cross(r_GE(:,1),vv);
        ee = cross(vv,hh)./mu - r_GE(:,1)./norm(r_GE(:,1));
        ev = ee/norm(ee).*(a*(1-norm(ee)));
        plot3(linspace(0,ev(1),100),linspace(0,ev(2),100),linspace(0,ev(3),100),'y','LineWidth',2)
        plot3(ev(1),ev(2),ev(3),'oy','LineWidth', 6)

    end
end
grid on
hold on

axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')