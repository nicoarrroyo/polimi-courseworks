%% configure paths
script_path = fileparts(mfilename("fullpath"));
backs = strfind(script_path, "\"); labs_dir = script_path(1:backs(end));
addpath([labs_dir '\student_functions']); addpath([labs_dir '\lib']);
clear; close all; clc;

%% constants / initial conditions
w_E = deg2rad(15.04) / 3600; % earth rotation velocity [rad s^-1]
mu_E = astroConstants(13); % earth gravitational parameter [km^3 s^-2]
earth_img = imread("EarthTexture.jpg");

%% %%% generic orbit case (k = 12, m = 1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kep1 = [
    8350, ...               % a: semi-major axis
    0.1976, ...             % e: eccentricity
    60, ...                 % i: inclination
    270, ...                % Omega: right ascension of the ascending node
    45, ...                 % omega: argument of perigee
    230];                   % TA: true anomaly]
[r01, v01] = kep2car( ...
    kep1(1), ...            % a: semi major axis (repeating)
    kep1(2), ...            % e: eccentricity
    kep1(3), ...            % i: inclination
    kep1(4), ...            % Omega: right ascension of the ascending node
    kep1(5), ...            % omega: argument of perigee
    kep1(6));               % TA: true anomaly

k1 = 12; m1 = 1;
a1r = (mu_E / ((w_E * k1 / m1) ^ 2)) ^ (1/3); % repeating a
[r01r, v01r] = kep2car( ...
    a1r, ...                % a: semi major axis (repeating)
    kep1(2), ...            % e: eccentricity
    kep1(3), ...            % i: inclination
    kep1(4), ...            % Omega: right ascension of the ascending node
    kep1(5), ...            % omega: argument of perigee
    kep1(6));               % TA: true anomaly

orbits = 16; % number of orbits to propogate for
title = "Generic Orbit";
lab2_ground_track( r01, v01, orbits, title )
title = "Generic Orbit (Repeating)";
lab2_ground_track( r01r, v01r, orbits, title )
disp("generic orbit complete")

%% kep2car car2kep conversion tests
kep_pre = [
    8350, ...               % a: semi-major axis
    0.1976, ...             % e: eccentricity
    60, ...                 % i: inclination
    270, ...                % Omega: right ascension of the ascending node
    45, ...                 % omega: argument of perigee
    230];                   % TA: true anomaly]
[r_test, v_test] = kep2car( ...
    kep_pre(1), ...         % a: semi major axis (repeating)
    kep_pre(2), ...         % e: eccentricity
    kep_pre(3), ...         % i: inclination
    kep_pre(4), ...         % Omega: right ascension of the ascending node
    kep_pre(5), ...         % omega: argument of perigee
    kep_pre(6));            % TA: true anomaly
[...
    a_post, ...
    e_post, ...
    i_post, ...
    Omega_post, ...
    omega_post, ...
    TA_post] = car2kep(r_test, v_test);
kep_post = [...
    a_post, ...
    e_post, ...
    i_post, ...
    Omega_post, ...
    omega_post, ...
    TA_post];
diff = kep_pre - kep_post;
disp("total error due to kep-car-kep conversion")
disp(diff)