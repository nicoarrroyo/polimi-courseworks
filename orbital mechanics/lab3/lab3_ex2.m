%% configure paths
script_path = fileparts(mfilename("fullpath"));
backs = strfind(script_path, "\"); labs_dir = script_path(1:backs(end));
addpath([labs_dir '\student_functions']); addpath([labs_dir '\lib']);
clear; close all; clc;

% Exercise 2: Orbit Transfer Problem
%% 1. Compute the initial and final states in Cartesian coordinates
% initial state
[r1, v1] = kep2car( ... % initial state vector [km, km s^-1]
    12500, ...          % a: semi major axis (repeating)
    0, ...              % e: eccentricity
    0, ...              % i: inclination
    0, ...              % Omega: right ascension of the ascending node
    0, ...              % omega: argument of perigee
    120);               % TA: true anomaly
r1 = r1'; v1 = v1'; % transpose to match Lambert function dimensions

% final state
[r2, v2] = kep2car( ... % final state vector [km, km s^-1]
    9500, ...           % a: semi major axis (repeating)
    0.3, ...            % e: eccentricity
    0, ...              % i: inclination
    0, ...              % Omega: right ascension of the ascending node
    0, ...              % omega: argument of perigee
    250);               % TA: true anomaly
r2 = r2'; v2 = v2'; % transpose to match Lambert function dimensions

%% 2. Solve Lambert's problem for the transfer arc
tof = 3300; % time of flight [s]
t1 = 0; % initial time [s]
mu_E = astroConstants(13); % earth gravitational parameter [quello che e`]

[A, P, E, ERROR, vt1, vt2, TPAR, THETA] = ...
    lambertMR( r1, r2, tof, mu_E, 0, 0, 0, 0 );

%% 3. Compute the total cost of the manoeuvre (Δv1 + Δv2)
dv1 = norm(abs(vt1 - v1));
dv2 = norm(abs(vt2 - v2));
dvtot = dv1 + dv2;

%% 4. Propagate the transfer arc from t1 to t2
yt = [r1 vt1]; % transfer arc initial state vector

% time span
t2 = t1 + tof;
tspan = linspace(t1, t2, 5000);

% set ODE solver conditions
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% integrate
[T, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, yt, options);
r = Y(:, 1:3);
v = Y(:, 4:6);

%% 5. Plot the initial and final orbits, and the transfer arc
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