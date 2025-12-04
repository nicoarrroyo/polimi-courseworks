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

%% %%% molniya orbit case (k = 2, m = 1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kep2 = [
    26600, ...               % a: semi major axis
    0.74, ...               % e: eccentricity
    63.4, ...               % i: inclination
    50, ...                 % Omega: right ascension of the ascending node
    280, ...                % omega: argument of perigee
    0];                     % TA: true anomaly
[r02, v02] = kep2car( ...
    kep2(1), ...            % a: semi major axis (repeating)
    kep2(2), ...            % e: eccentricity
    kep2(3), ...            % i: inclination
    kep2(4), ...            % Omega: right ascension of the ascending node
    kep2(5), ...            % omega: argument of perigee
    kep2(6));               % TA: true anomaly

k2 = 2; m2 = 1;
a2r = (mu_E / ((w_E * k2 / m2) ^ 2)) ^ (1/3); % a from km
[r02r, v02r] = kep2car( ...
    a2r, ...                % a: semi major axis (repeating)
    kep2(2), ...            % e: eccentricity
    kep2(3), ...            % i: inclination
    kep2(4), ...            % Omega: right ascension of the ascending node
    kep2(5), ...            % omega: argument of perigee
    kep2(6));               % TA: true anomaly

orbits = 30; % number of orbits to propogate for
title = "Molniya Orbit";
lab2_ground_track( r02, v02, orbits, title )
title = "Molniya Orbit (Repeating)";
lab2_ground_track( r02r, v02r, orbits, title )
disp("molniya orbit complete")

%% %%% three circular LEO orbits case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kep3 = [
    [7171.010, 7171.010, 7171.010];     % a: semi major axis
    [0, 0, 0];                          % e: eccentricity
    [0, 30, 98];                        % i: inclination
    [0, 0, 0];                          % Omega: right ascension of the ascending node
    [40, 40, 40];                       % omega: argument of perigee
    [0, 0, 0];];                        % TA: true anomaly
r03 = zeros(3, 3); v03 = zeros(3, 3);
for i = 1:3
    [r03(i, :), v03(i, :)] = kep2car( ...
        kep3(1, i), ...            % a: semi major axis (repeating)
        kep3(2, i), ...            % e: eccentricity
        kep3(3, i), ...            % i: inclination
        kep3(4, i), ...            % Omega: right ascension of the ascending node
        kep3(5, i), ...            % omega: argument of perigee
        kep3(6, i));               % TA: true anomaly
end

k3 = [20, 29, 15]; m3 = [2, 2, 1];
a3r = (mu_E ./ ((w_E .* k3 ./ m3) .^ 2)) .^ (1/3); % a from km
r03r = zeros(3, 3); v03r = zeros(3, 3);
for i = 1:3
    [r03r(i, :), v03r(i, :)] = kep2car( ...
        a3r(i), ...                % a: semi major axis (repeating)
        kep3(2, i), ...            % e: eccentricity
        kep3(3, i), ...            % i: inclination
        kep3(4, i), ...            % Omega: right ascension of the ascending node
        kep3(5, i), ...            % omega: argument of perigee
        kep3(6, i));               % TA: true anomaly
end

for i = 1:3
    title = "Circular LEO " + num2str(i);
    lab2_ground_track( r03(i, :), v03(i, :), orbits, title )

    title = "Circular LEO " + num2str(i) + " (Repeating)";
    lab2_ground_track( r03r(i, :), v03r(i, :), orbits, title )
    disp("circular LEO orbit " + i + " complete")
end

%% %%% geostationary and geosynchronous orbit case %%%%%%%%%%%%%%%%%%%%%%%%
[r04, v04] = kep2car( ...
    42166.167, ...
    0, ...
    0, ...
    0, ...
    0, ...
    20);
orbits = 50;
title = "Geostationary Orbit";
lab2_ground_track( r04, v04, orbits, title )
disp("GEO ground track complete")

[r05, v05] = kep2car( ...
    42166.167, ...
    0, ...
    25, ...
    0, ...
    0, ...
    20);
orbits = 10;
title = "Geosynchronous Orbit";
lab2_ground_track( r05, v05, orbits, title )
disp("GSO ground track complete")

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
disp(diff)