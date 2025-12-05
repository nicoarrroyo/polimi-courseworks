%% configure paths
script_path = fileparts(mfilename("fullpath"));
backs = strfind(script_path, "\"); labs_dir = script_path(1:backs(end));
addpath([labs_dir '\student_functions']); addpath([labs_dir '\lib']);
clear; close all; clc;

% Exercise 1: State reconstruction problem
%% 1. Write a script to solve Lambert's problem for the values of r(t1), 
% r(t2) and time of flight tof = t2 - t1 given:

r1 = [-21800 37900 0]; % initial position [km]
r2 = [27300 27700 0]; % final position [km]
tof = (15 * 3600) + (6 * 60) + 40; % time of flight [s]
mu_E = astroConstants(13);

[A, P, E, ERROR, v1, v2, TPAR, THETA] = ...
    lambertMR( r1, r2, tof, mu_E, 0, 0, 0, 0 );

%% 2. Propagate and plot the resulting orbit
y1 = [r1 v1];

% time span
a = 1/(2/norm(r1) - dot(v1,v1)/mu_E); % [km]
T = 2*pi*sqrt(a^3/mu_E); % [s]
tspan = linspace(0, T, 5000);

% set ODE solver conditions
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% integrate
[T, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y1, options);
r = Y(:, 1:3);
v = Y(:, 4:6);

% %%% a. plot the orbit over 1 period T %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get texture of earth
earth_img = imread("EarthTexture.jpg");
[x_earth, y_earth, z_earth] = sphere(50);
R_e = astroConstants(23); % earth radius [km]
x_earth = R_e * x_earth;
y_earth = R_e * y_earth;
z_earth = -R_e * z_earth;

% plot
figure("Name", "Orbit Plot");
plot3(r(:, 1), r(:, 2), r(:, 3));
hold on
surface(x_earth, y_earth, z_earth, "FaceColor", "texturemap", ...
    "CData", earth_img, "EdgeColor", "none")
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Two-body problem orbit");
scatter3(r1(1), r1(2), r1(3), "filled", "MarkerFaceColor", "r");
scatter3(r2(1), r2(2), r2(3), "filled", "MarkerFaceColor", "g");
legend("orbit", "earth", "initial point", "final point");
axis equal; grid on;
hold off