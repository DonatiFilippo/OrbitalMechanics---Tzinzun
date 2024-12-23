function plotOrbit(a, e, i, OM, om, theta, point, th1, th2, tratteggio)

% Function to plot a 3D orbit in the inertial frame
%
%
% plotorbit(a, e, i, OM, om, theta, point, th1, th2, tratteggio)
%
% Input arguments
% -------------------------------------------------------
% a     [1x1]     semimajor axis        [km]
% e     [1x1]     eccentricity          [-]
% i     [1x1]     inclination           [rad]
% OM    [1x1]     RAAN                  [rad]
% om    [1x1]     pericenter anomaly    [rad]
% theta [1x1]     true anomaly of the labled point  [rad]
% point [string]  lable of the pint     
% th1   [1x1]     true anomaly of the start of the curve line [rad]
% th2   [1x1]     true anomaly of the end of the curve line [rad]
% tratteggio [boolean] Set to "si" to have dashed line 
% -------------------------------------------------------

    if nargin == 7 || nargin == 5
        th1=0;
        th2=2*pi;
        tratteggio = 0;
    end

    dth=1e-2;

    % Plot the Orbit: 
    th2 = th2 + 2*pi*(th2 < th1);
    cax=newplot();

    % Dashed Orbit Plot:
    if tratteggio == 'si'
        rvet=[];
        for th=0:dth:2*pi
            [rr, ~] = parOrb2rv(a, e, i, OM, om, th);
            rvet = [rvet rr];
        end
    
        plot3(cax, rvet(1,:), rvet(2,:), rvet(3,:),'k--','LineWidth',0.1);
        quiver3(0,0,0,2*rvet(1,1),2*rvet(2,1),2*rvet(3,1),'LineWidth',1);
    end

    % Plot of the continuos line between th1 and th2
    rvet=[];
    for th=th1:dth:th2
        [rr, ~] = parOrb2rv(a, e, i, OM, om, th);
        rvet = [rvet rr];
    end

    plot3(cax, rvet(1,:), rvet(2,:), rvet(3,:),'LineWidth',1.7);
    grid on;
    axis equal;


    if nargin > 5
        [rr, ~] = parOrb2rv(a, e, i, OM, om, theta);
    
        plot3(cax, rr(1), rr(2), rr(3),'ok','MarkerFaceColor','k');
    
        text(rr(1), rr(2), rr(3), point, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 15);
    end
end