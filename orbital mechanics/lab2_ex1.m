clear; close all; clc;

%% constants / initial conditions
w_E = 15.04; % earth rotation velocity [deg hr^-1]
mu_E = astroConstants(13); % earth gravitational parameter [km^3 s^-2]
R_e = astroConstants(23); % earth radius [km]
earth_img = imread("EarthTexture.jpg");

%% %%% generic orbit case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% intial conditions
% e = 0.1976; % eccentricity [-]
i = 60; % inclination [deg]
% omega = 270; % [deg]
% w = 45; % [deg]
% f0 = 230; % [deg]
orbits = 3.25; % number of orbits to propogate for [-]

r0 = [-4578.219 -801.084 -7929.708]; % position vector [km]
v0 = [0.800 -6.037 1.385]; % velocity vector [km s^-1]
y0 = [r0 v0]; % state vector

a = 1 / (2 / norm(r0) - dot(v0,v0) / mu_E); % semi-major axis [km]

%% time span
period = 2 * pi * sqrt(a^3 / mu_E); % orbital period [s]
tspan = linspace(0, orbits*period, 10000);

%% orbit propogation
% set solver conditions
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% integrate
[T, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y0, options);
r = Y(:, 1:3); % FOR TROUBLESHOOTING

%% ground track calculation and plotting
% calculation
[long, lat] = ground_track(Y, T, w_E); % longitude and latitude [radians]

% plotting
figure("Name", "Generic Orbit Ground Track Plot")
image([-180 180], [90 -90], earth_img);

hold on
plot(rad2deg(long), rad2deg(lat), ".")
xline(rad2deg(long(1)), "r");
xline(rad2deg(long(end)), "g");

title("ground track plot")
xlabel("longitude [deg]"); ylabel("latitude [deg]");
grid on;
ylim([-90 90]); xlim([-180 180]);
legend("ground track", "initial long", "final long")
set(gca, "Ydir", "normal")
hold off

%% %%% molniya orbit case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%% three circular LEO orbits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

