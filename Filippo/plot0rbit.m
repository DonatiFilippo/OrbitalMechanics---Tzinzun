function plot0rbit(a, e, i, OM, om, th0, thf, dth, mu, mod)
% ORBIT REPRESENTATION
%
%--------------------------------------------------------------------------
%INPUT ARGUMENTS:
%a        [1x1]  semi-major axis                     [km]
%e        [1x1]  eccentricity                        [-]
%om       [1x1]  argument of periapsis               [rad]
%OM       [1x1]  longitude of ascending node         [rad]
%th0      [1x1]  starting true anomaly               [rad]
%thf      [1x1]  final true anomaly                  [rad]
%dth      [1x1]  true anomaly resolution             [rad]
%mu       [1x1]  gravitational parameter             [km^3/s^2]
%mod      [1x1]  discrimination parameter            [-]
%
%--------------------------------------------------------------------------
%OUTPUT:
%Graphic representation of the orbit from th0 to thf in 3D.
%--------------------------------------------------------------------------
%AUTHORS: 
% Eleonora Domenichelli, Elena M. Chianese, Riccardo Coppola
%
%--------------------------------------------------------------------------


r_GE=[];

for t=th0:dth:thf
[rr] = parorb2rv (a, e, i, OM, om, t, mu);
r_GE=[r_GE rr];
end

if nargin == 9
plot3(r_GE(1,:), r_GE(2,:), r_GE(3,:),'LineWidth',2)

else
    if mod == 1 %dashed orbit in black
        plot3(r_GE(1,:), r_GE(2,:), r_GE(3,:),'k--')

    elseif mod == 2 %colored continuos line orbit with eccentricity vector and pericenter
        plot3(r_GE(1,:), r_GE(2,:), r_GE(3,:),'LineWidth',2)

        [r_GE(:,1),vv] = parorb2rv (a, e, i, OM, om,th0, mu);
        hh = cross(r_GE(:,1),vv);
        ee = cross(vv,hh)./mu - r_GE(:,1)./norm(r_GE(:,1));
        ev = ee/norm(ee).*(a*(1-norm(ee)));
        plot3(linspace(0,ev(1),100),linspace(0,ev(2),100),linspace(0,ev(3),100),'y','LineWidth',2)
        plot3(ev(1),ev(2),ev(3),'oy','LineWidth', 6)

    elseif mod == 3 %dashed orbit in black with eccentricity vector and pericenter
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