%% Main script for Planetary Explorer Mission assignment %%
%
% Group number : 2427
%
%--------------------------------------------------------------------------
% Created and maintained by : 
%
% Azevedo Da Silva Esteban
% Gavidia Pantoja Maria Paulina
% Donati Filippo 
% Domenichelli Eleonora
% 
%--------------------------------------------------------------------------
% LAST UPDATE: 22-12-2024
%

% 2427 / Earth / J2 & Moon / STARTING DATE OF PROPAGATION? 

%% SCRIPT INITIALISATION
clc
clear variables
close all

%% PRESENTATION
fprintf('----------------------------------------------\n');
fprintf('   PLANETARY EXPLORER MISSION ASSIGNMENT      \n');
fprintf('----------------------------------------------\n');
fprintf('\n');

%% Important Constants
mu_E = astroConstants(13);
mu_Moon = astroConstants(20); 
R_E = astroConstants(23);
J2 = astroConstants(9);
w_E = deg2rad(15.04) * (1/3600);

%% ORBIT CHARACTERISATION
a = 23300; % Semi-major axis [km]
e = 0.523; % Eccentricity [-]
i = deg2rad(58.65); % Inclination [rad]
OM = 0; % RAAN [rad]
om = 0; % Argument of pericenter [rad]
theta0 = 0; % True Anomaly [rad]

date = [2024, 12, 22, 23, 20, 30];

Earth_3D;
hold on
% plotOrbit(a, e, i, OM, om);

k = 13;
m = 3;

%% Ground Track

N = 100000;

T = 2*pi*sqrt(a^3/mu_E);
t1 = linspace(0, T, N);
t2 = linspace(0, 24*60*60, N);
t3 = linspace(0, 30*T, N);

[rr, vv] = kep2car(a, e, i, OM, om, theta0, mu_E);
y0g = [a; e; i; OM; om; theta0];
y0 = [rr; vv];
th = thetaG(2024, 12, 22, 23, 20, 30,w_E);
groundTrack2(y0,t1,th,mu_E,w_E);
groundTrack2(y0,t2,th,mu_E,w_E);
groundTrack2(y0,t3,th,mu_E,w_E);

n = w_E*(k/m);
a_rep = power(mu_E/(n^2),1/3);
Tr = 2*pi*sqrt(a_rep^3/mu_E);
t1r = linspace(0, Tr, N);
t2r = linspace(0, 3*24*60*60, N);
[rrr, vvr] = kep2car(a_rep, e, i, OM, om, theta0, mu_E);
y0r = [rrr; vvr];
groundTrack2(y0r,t1r,th,mu_E,w_E);
groundTrack2(y0r,t2r,th,mu_E,w_E);
groundTrack2(y0r,t3,th,mu_E,w_E);

st = date2mjd2000(date);


fprintf('ne mancano solo 6\n');
fprintf('\n');

groundTrackJ2(y0,t1,th,mu_E,w_E, J2, R_E, mu_Moon, st);
groundTrackJ2(y0,t2,th,mu_E,w_E, J2, R_E, mu_Moon, st);
groundTrackJ2(y0,t3,th,mu_E,w_E, J2, R_E, mu_Moon, st);

fprintf('quasi arrivato\n');
fprintf('\n');

groundTrackJ2(y0r,t1r,th,mu_E,w_E, J2, R_E, mu_Moon, st);
groundTrackJ2(y0r,t2r,th,mu_E,w_E, J2, R_E, mu_Moon, st);
groundTrackJ2(y0r,t3,th,mu_E,w_E, J2, R_E, mu_Moon, st);

%% Le mucche
t3m = linspace(0, 30*T, N);
fprintf('Le mucche fanno mu, ma una fa mu Moon\n');
groundTrackJ2G(y0g,t1,th,mu_E,w_E, J2, R_E, mu_Moon, st);
groundTrackJ2G(y0g,t2,th,mu_E,w_E, J2, R_E, mu_Moon, st);
groundTrackJ2G(y0g,t3m,th,mu_E,w_E, J2, R_E, mu_Moon, st);

%% Punto 4
tv = linspace(0, 100*T, N);
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[tu, yu] = ode113(@(t, y) gaussFun(t, y, accPert(t, y,mu_E, J2, R_E, mu_Moon,st),mu_E), tv, y0g, options);

l = size(tv,2);
error = zeros(3,l);
for i = 1:l
    rrv = kep2car(yu(i, 1), yu(i, 2), yu(i, 3), yu(i, 4), yu(i, 5), yu(i, 6), mu_E);
    aNora(:, i) = a3B_M(tu(i),rrv,mu_Moon, date');
    afil(:, i) = accMoon(tu(i),yu(i,:)', mu_Moon, mu_E, st);
end 


% figure
% plot(tu,yu(:, 1));
% title("a");
% figure
% plot(tu,yu(:, 2));
% title("e");
% figure
% plot(tu,yu(:, 3));
% title("i");
% figure
% plot(tu,yu(:, 4));
% title("OM");
% figure
% plot(tu,yu(:, 5));
% title("om");
% figure
% plot(tu,yu(:, 6));
% title("theta");
%% 
tv = linspace(0, 50*T, N);
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[tuc, keple] = ode113(@(t, y) twoBodyProblemPert(t, y, mu_E, J2, R_E, mu_Moon, st), tv, y0, options);
l = size(tuc,1);
s = zeros(l, 6);
for i = 1:l
    [s(i,1),s(i,2), s(i,3), s(i,4), s(i,5), s(i,6)] = rv2parorb(keple(i, 1:3), keple(i, 4:6), mu_E);
end
s(:,4) = unwrap(s(:,4));
s(:,5) = unwrap(s(:,5));
s(:,6) = unwrap(s(:,6));

figure
plot(s(:, 1));
title("a_err");
figure
plot(s(:, 2));
title("e_err");
figure
plot(s(:, 3));
title("i_err");
figure
plot(s(:, 4));
title("OM_err");
figure
plot(s(:, 5));
title("om_err");
figure
plot(s(:, 6));
title("theta_err");
fprintf('Tnh\n');
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[tuc, yut] = ode113(@(t, y) gaussTNH(t, y, accPertTNH(t, y, mu_E, J2, R_E, mu_Moon, st),mu_E), tv, y0g, options);

%% 
figure

plot(tuc/T,yut(:, 1));
grid on
title("a");
figure
plot(tuc/T,yut(:, 2));
title("e");
grid on
figure

plot(tuc/T,yut(:, 3));
title("i");
grid on
figure

plot(tuc/T,rad2deg(yut(:, 4)));
title("OM");
grid on
figure

plot(tuc/T,yut(:, 5));
grid on
title("om");
figure

plot(tuc/T,yut(:, 6));
grid on
title("theta");

%%
Earth_3D;
hold on
scatter3( keple(:,1), keple(:,2), keple(:,3), 10, tuc./T, 'filled');
colormap;
a = colorbar;
a.Title.String = "Orbital Periods [T]";
clim([0 50]);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;

%% Filtering (1000T)
as = s(:,1);
filterA = movmean(as, N/500);
figure
plot(tv/T, filterA)

es = s(:, 2);
filterE = movmean(es, N*2/33);
figure
plot(tv/T, filterE)



%% Filter I
figure
is = s(:, 3);
plot(tv/T, is, 'LineWidth',0.2)
hold on


filterI = movmean(is, N/500);
plot(tv/T, filterI, 'Color', [1, 0, 0], 'LineWidth', 0.6)

filterIs = movmean(is, N*2/33);
plot(tv/T, filterIs, 'Color', [0, 1, 0], 'LineWidth', 1)
title('Inclination')
%%
OMs = s(:, 4);
filterOM = movmean(OMs, [0, N/500]);
figure
plot(tv/T, filterOM)

oms = s(:, 5);
filterom = movmean(oms, N/500);
figure
plot(tv/T, filterom)

ths = s(:, 6);
filterTH = movmean(ths, N/500);
figure
plot(tv/T, filterTH)


