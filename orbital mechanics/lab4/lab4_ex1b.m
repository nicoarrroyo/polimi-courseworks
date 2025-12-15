%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
labs_d = cd(1:backs(end)); addpath([labs_d '\student_functions']); 
addpath([labs_d '\lib']); addpath([labs_d '\lib' '\timeConversion']);
clear; close all; clc;

% Exercise 1b
% Design a fly-by around the Earth for a fixed location of the incoming 
% asymptote and different impact parameters.
%% 1. choose a location for the incoming asymptote
% --- constants ---
mu_E                = astroConstants(13);   % earth grav. param. [km^3 s^-2]
mu_sun              = astroConstants(4);    % sun grav. param. [km^3 s^-2]
AU                  = astroConstants(2);    % astronomical unit [km]
R_E                 = AU .* [1; 0; 0;];     % earth orbital radius [km]
r_E                 = astroConstants(23);   % earth mean radius [km]

% --- decided terms ---
impact_param_mags = [8200; 10200; 12200; 14200; 16200;];
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);
steps = 1500;
tspan = linspace(0, 10000, steps);

% --- initial conditions ---
v_inf_minus = [15.1; 0; 0;]; % velocity before fly-by [km s^-1]

% --- calculated terms ---
a          = -mu_E / (norm(v_inf_minus) ^ 2);   % semi-major axis [km]
turn_angle = 2 .* atan(-a ./ impact_param_mags); % turn angle, delta [rad]
ecc        = 1 / sin(turn_angle / 2);           % eccentricity [-]
r_p        = a * (1 - ecc);                     % pericentre radius [km]

%% 2. solve and plot the 2d hyperbola for all values of impact parameter
% --- solve the 2d hyperbola (geocentric) for all values of impact parameter ---
dv = zeros(length(impact_param_mags), 3);
r = zeros(length(impact_param_mags), steps, 3);

% --- create figure ---
figure("Name", "Different Impact Parameters"); hold on;

% --- start orbit propagation loop logic loop ---
for i = 1:length(impact_param_mags)
    impact_param = [0; impact_param_mags(i); 0;];

    % --- build inertial frame ---
    u = cross(-impact_param, v_inf_minus); % vector normal to plane of hyperbola
    u_hat = u / norm(u); % unit vector normal to plane of hyperbola
    v_inf_plus = rodrigues(v_inf_minus, u_hat, turn_angle(i));
    dv(i, :) = v_inf_plus - v_inf_minus;
    
    % --- simulate fly-by with 2bp dynamics ---
    y = [[-2*r_E; impact_param_mags(i); 0;] v_inf_minus];
    [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y, options);
    r(i, :, :) = Y(:, 1:3);

    % --- plot all impact parameter fly-by trajectories ---
    plot3(r(i, :, 1), r(i, :, 2), r(i, :, 3));
end

% --- plot earth surface texture ---
earth_img = imread("EarthTexture.jpg");
[x_earth, y_earth, z_earth] = sphere(50);
x_earth = r_E * x_earth; y_earth = r_E * y_earth; z_earth = -r_E * z_earth;
surface(x_earth, y_earth, z_earth, "FaceColor", "texturemap", ...
    "CData", earth_img, "EdgeColor", "none")

% --- finish up plot properties ---
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Earth Fly-by Trajectories for a Range of Impact Parameters");
grid on; axis equal;
xlim([-2*r_E inf]);
hold off;

% --- output final results ---
% disp("=== LAB 3 EX. 1b RESULTS ===")
% fprintf("\n--- COMMON RESULTS ---\n")
% fprintf("       turn angle: %.4f deg\n", rad2deg(turn_angle));
% fprintf("pericentre radius: %.4f km s^-1\n", r_p);
% fprintf("   fly-by delta-v: %.4f km s^-1\n", dv);
% fprintf("  semi-major axis: %.4f km\n", a);
% fprintf("     eccentricity: %.4f km\n", ecc);

%% 3. compute V_minus and the incoming heliocentric arc for each Δ
% --- compute V_minus ---
V_E = sqrt(mu_sun / norm(R_E)) * [0; 1; 0;]; % earth heliocentric velocity [km s^-1]
V_minus = V_E + v_inf_minus;

% --- propagate the incoming heliocentric arc ---
t_span_helio = linspace(0, 2.628 * 10^6, steps); % 1 month in seconds

disp("propagating incoming heliocentric arc")
for i = 1:length(impact_param_mags)
    y = [[AU; 0; 0;] V_minus];
    [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), -t_span_helio, y, options);
end
R_minus = Y(:, 1:3);
disp("complete!")

%% 4. compute v_inf_plus, V_plus, and the outgoing heliocentric arc for each Δ
% --- compute V_plus ---
v_inf_plus = zeros(length(impact_param_mags), 3);
V_plus = zeros(length(impact_param_mags), 3);
for i = 1:length(impact_param_mags)
    v_inf_plus(i, :) = v_inf_minus + dv(i, :)';
    V_plus(i, :) = V_E + v_inf_plus(i, :)';
end

% --- propagate the outgoing heliocentric arc ---
R_plus = zeros(length(impact_param_mags), steps, 3);

disp("propagating outgoing heliocentric arc")
for i = 1:length(impact_param_mags)
    y = [R_E V_plus(i, :)'];
    [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), 23*t_span_helio, y, options);
    R_plus(i, :, :) = Y(:, 1:3);
end
disp("complete!")

%% 5. plot the heliocentric trajectory
% --- create figure ---
figure("Name", "Heliocentric Different Impact Parameters"); hold on;

% --- plot position of earth during fly-by ---
% also propagate earth position before and after fly-by
y = [R_E V_E];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), -t_span_helio, y, options);
R_E_minus = Y(:, 1:3);
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), 11*t_span_helio, y, options);
R_E_plus = Y(:, 1:3);
plot3(R_E_minus(:, 1), R_E_minus(:, 2), R_E_minus(:, 3), "g--");
plot3(R_E_plus(:, 1), R_E_plus(:, 2), R_E_plus(:, 3), "g");
scatter3(R_E(1), R_E(2), R_E(3), "filled", "g");

% --- plot heliocentric trajectories for each impact parameter ---
plot3(R_minus(:, 1), R_minus(:, 2), R_minus(:, 3), "b");
for i = 1:1:length(impact_param_mags)
    plot3(R_plus(i, :, 1), R_plus(i, :, 2), R_plus(i, :, 3));
end

% --- plot sun ---
scatter3(0, 0, 0, "filled", "y");

% --- finish up plot properties ---
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Earth Fly-by Trajectories for a Range of Impact Parameters");
grid on; axis equal;
legend( ...
    "", "Earth Trajectory", "", ...
    "Incoming Fly-By Trajectory", ...
    "Δ=8200 Outgoing Fly-By Trajectory", ...
    "Δ=10200", "Δ=12200", "Δ=14200", "Δ=16200", "Sun")
hold off;
