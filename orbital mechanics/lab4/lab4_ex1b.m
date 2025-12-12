%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
labs_d = cd(1:backs(end)); addpath([labs_d '\student_functions']); 
addpath([labs_d '\lib']); addpath([labs_d '\lib' '\timeConversion']);
clear; close all; clc;

% Exercise 1a
% Design a fly-by around the Earth for a fixed location of the incoming 
% asymptote and different impact parameters.
%% 1. choose a location for the incoming asymptote
% --- constants ---
mu_E                = astroConstants(13);   % earth grav. param. [km^3 s^-2]
mu_sun              = astroConstants(4);    % sun grav. param. [km^3 s^-2]
AU                  = astroConstants(2);    % astronomical unit [km]
R_E                 = AU .* [0; 1; 0;];     % earth orbital radius [km]
r_E                 = astroConstants(23);   % earth mean radius [km]

% --- decided terms ---
impact_param_mags = [8200; 10200;];%; 12200; 14200; 16200; 18200;];
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);
steps = 500;
tspan = linspace(0, 10000, steps);

% --- initial conditions ---
v_inf_minus = [15.1; 0; 0;]; % velocity before fly-by [km s^-1]

% --- calculated terms ---
a           = - mu_E / (norm(v_inf_minus) ^ 2); % semi-major axis [km]
turn_angle  = 2 .* atan(- a ./ impact_param_mags); % turn angle, delta [rad]
ecc         = 1 / sin(turn_angle / 2);          % eccentricity [-]
r_p         = a * (1 - ecc);                    % pericentre radius [km]

%% 2. solve and plot the 2d hyperbola for all values of impact parameter
% --- solve the 2d hyperbola for all values of impact parameter ---
dv = zeros(length(impact_param_mags), 3);
r = zeros(length(impact_param_mags), steps, 3);
for i = 1:length(impact_param_mags)
    impact_param = [0; impact_param_mags(i); 0;];

    u = cross(-impact_param, v_inf_minus); % vector normal to plane of hyperbola
    u_hat = u / norm(u); % unit vector normal to plane of hyperbola
    v_inf_plus = rodrigues(v_inf_minus, u_hat, turn_angle(i));
    dv(i, :) = v_inf_plus - v_inf_minus;
    
    y = [[-2*r_E; impact_param_mags(i); 0;] v_inf_minus];
    [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y, options);
    r(i, :, :) = Y(:, 1:3);
end

% --- create figure ---
figure("Name", "Different Impact Parameters"); hold on;

% --- plot all impact parameter fly-by trajectories
plot3(r(:, :, 1), r(:, :, 2), r(:, :, 3));

% --- plot earth surface texture ---
earth_img = imread("EarthTexture.jpg");
[x_earth, y_earth, z_earth] = sphere(50);
x_earth = r_E * x_earth; y_earth = r_E * y_earth; z_earth = -r_E * z_earth;
surface(x_earth, y_earth, z_earth, "FaceColor", "texturemap", ...
    "CData", earth_img, "EdgeColor", "none")

% --- finish up plot properties ---
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Earth Fly-by Trajectories for a Range of Impact Parameters");

axis equal; grid on;
xlim([-2*r_E 16*r_E]);
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
V_E = sqrt(mu_sun / R_E)'; % earth heliocentric velocity [km s^-1]
V_minus = V_E + v_inf_minus;

%% 4. compute v_inf_plus, V_plus, and the outgoing heliocentric arc for each Δ
v_inf_plus = v_inf_minus + dv';
V_plus = V_E + v_inf_plus;

%% 5. plot the heliocentric trajectory
% --- create figure ---