%% configure paths
script_path = fileparts(mfilename("fullpath"));
backs = strfind(script_path, "\"); labs_dir = script_path(1:backs(end));
addpath([labs_dir '\student_functions']); addpath([labs_dir '\lib']);
clear; close all; clc;

%% constants / initial conditions
w_E = deg2rad(15.04) / 3600; % earth rotation velocity [rad s^-1]
mu_E = astroConstants(13); % earth gravitational parameter [km^3 s^-2]
J2 = astroConstants(9); % second zonal harmonic
R_e = astroConstants(23); % earth radius [km]
earth_img = imread("EarthTexture.jpg");
a = 11376; e = 0.4178; i = 66; Omega = 40; omega = 20; TA = 23;

%% %%% question 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kep = [...
    a, ...
    e, ...
    i, ...
    Omega, ...
    omega, ...
    TA];
[r0, v0] = kep2car( ...
    kep(1), ...             % a: semi major axis (repeating)
    kep(2), ...             % e: eccentricity
    kep(3), ...             % i: inclination
    kep(4), ...             % Omega: right ascension of the ascending node
    kep(5), ...             % omega: argument of perigee
    kep(6));                % TA: true anomaly

fprintf("===== QUESTION 1 RESPONSE =====\n")
disp("INITIAL POSITION VECTOR:")
fprintf("%.2f\n", r0)
disp("INITIAL VELOCITY VECTOR:")
fprintf("%.2f\n", v0)

%% %%% question 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up ODE solver conditions
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);
tspan = linspace(0, 3*24*3600, 1000); % time span for integration
y0 = [r0; v0]; % initial state vector

% integrate for unperturbed 2bp
[T_2bp, Y_2bp] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y0, options);
[long_2bp, lat_2bp] = ground_track(Y_2bp, T_2bp, w_E); % longitude and latitude [rad]

% integrate for j2-perturbed 2bp
[T_j2, Y_j2] = ode113(@(t,y) ode_2bp_j2(t,y,mu_E,J2,R_e), tspan, y0, options);
[long_j2, lat_j2] = ground_track(Y_j2, T_j2, w_E); % longitude and latitude [rad]

fprintf("\n===== QUESTION 2 RESPONSE =====\n")
% fprintf("2BP UNPERTURBED FINAL LONGITUDE: %.2e rad, %.2e deg\n", ...
%     long_2bp(end), rad2deg(long_2bp(end)));
% fprintf("2BP UNPERTURBED FINAL LATITUDE:  %.2e rad, %.2e deg\n", ...
%     lat_2bp(end), rad2deg(lat_2bp(end)));
% 
% fprintf("J2 PERTURBED FINAL LONGITUDE: %.2e rad, %.2e deg\n", ...
%     long_j2(end), rad2deg(long_j2(end)));
% fprintf("J2 PERTURBED FINAL LATITUDE:  %.2e rad, %.2e deg\n", ...
%     lat_j2(end), rad2deg(lat_j2(end)));

long_diff = abs(long_2bp(end) - long_j2(end));
lat_diff = abs(lat_2bp(end) - lat_j2(end));
fprintf("ABSOLUTE DIFFERENCE FINAL LONGITUDE: %.4f rad, %.2f deg\n", ...
    long_diff, rad2deg(long_diff));
fprintf("ABSOLUTE DIFFERENCE FINAL LATITUDE:  %.4f rad, %.2f deg\n", ...
    lat_diff, rad2deg(lat_diff));

%% %%% question 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 9; m = 3;
ar = (mu_E / ((w_E * k / m) ^ 2)) ^ (1/3); % repeating SMA
[rr, vr] = kep2car( ...
    ar, ...                 % a: semi major axis (repeating)
    kep(2), ...             % e: eccentricity
    kep(3), ...             % i: inclination
    kep(4), ...             % Omega: right ascension of the ascending node
    kep(5), ...             % omega: argument of perigee
    kep(6));                % TA: true anomaly

period = (2*pi) * sqrt(ar^3 / mu_E);
orbits = 3 * 24 * 3600 / period;
title = "Generic Orbit";

fprintf("\n===== QUESTION 3 RESPONSE =====\n")
lab2_ground_track( rr, vr, orbits, title )

%% %%% question 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ephemeris = [2460744.500000000, "2025-Mar-10 00:00:00.0000", ...
    3.361766760541229E+03, 1.003047961901085E+03, -6.023594759534308E+03, ...
    -6.457842573991658E+00, -6.913467679034817E-01, -3.744195427310673E+00];
r_3509 = [
    3.361766760541229E+03;
    1.003047961901085E+03;
    -6.023594759534308E+03;];
v_3509 = [
    -6.457842573991658E+00;
    -6.913467679034817E-01;
    -3.744195427310673E+00;];

fprintf("\n===== QUESTION 4 RESPONSE =====\n")
disp("SAT 3509 POSITION VECTOR:")
fprintf("%.2f\n", r_3509)
disp("SAT 3509 VELOCITY VECTOR:")
fprintf("%.2f\n", v_3509)