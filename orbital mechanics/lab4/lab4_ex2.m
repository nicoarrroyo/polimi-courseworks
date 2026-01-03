%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
labs_d = cd(1:backs(end)); addpath([labs_d '\student_functions']); 
addpath([labs_d '\lib']); addpath([labs_d '\lib' '\timeConversion']);
clear; close all; clc;

% Exercise 2
% Design a powered gravity assist around the Earth, given the incoming and
% outgoing heliocentric velocities, and the Earth's position
%% preamble
V_minus = [31.5; 5.2; 0;];    % incoming heliocentric velocity [km s^-1]
V_plus  = [36; 0; 0;];        % outgoing heliocentric velocity [km s^-1]
mu_E    = astroConstants(13); % earth gravitational parameter [km^2 s^-3]
mu_sun  = astroConstants(4);  % sun gravitational parameter [km^2 s^-3]
AU      = astroConstants(2);  % astronomical unit [km]
R_E     = [0; -1; 0;] * AU;   % earth orbit radius [km]
r_E     = astroConstants(23); % earth mean radius [km]

%% === part 1 ===
% Compute the velocities wrt the planet before and after the fly-by
V_E = sqrt(mu_sun / norm(R_E)) * [-R_E(2); R_E(1); 0;] / AU; % earth heliocentric velocity [km s^-1]
v_inf_minus = V_minus - V_E;
v_inf_plus = V_plus - V_E;

%% === part 2 ===
% Compute the turning angle
turn_angle = acos(dot(v_inf_minus, v_inf_plus) / (norm(v_inf_minus) * norm(v_inf_plus)));

%% === part 3 ===
% Solve the non-linear system for r_p and check its validity
h_atm = 500; % height of earth atmosphere from sea-level [km]

% --- solve the non-linear system for r_p ---
eq = @(r_p) turn_angle - ...
    asin(1 / (1 + (r_p * norm(v_inf_plus)^2) / mu_E)) - ...
    asin(1 / (1 + (r_p * norm(v_inf_minus)^2) / mu_E));
r_p = fzero(eq, r_E + h_atm);

% --- check its validity ---
r_p_crit = r_E + h_atm; % critical fly-by pericentre radius

if r_p < r_p_crit
    fprintf("\n!!! DANGER DANGER DANGER !!!\n")
    fprintf("\nfly-by pericentre radius might be too low!\n")
    fprintf("\nTERRAIN. PULL UP. TERRAIN. PULL UP.\n")
    fprintf("\n!!! DANGER DANGER DANGER !!!\n")
end

h_ga = r_p - r_E;

%% === part 4 ===
% Compute the velocities of the two hyperbolic arcs at pericentre and the 
% required delta-v @ pericentre
ecc_minus = 1 + (r_p * norm(v_inf_minus)^2 / mu_E);
a_minus = -mu_E / (norm(v_inf_minus)^2);
v_p_minus = sqrt(mu_E * (2/r_p - 1/a_minus));

ecc_plus = 1 + (r_p * norm(v_inf_plus)^2 / mu_E);
a_plus = -mu_E / (norm(v_inf_plus)^2);
v_p_plus = sqrt(mu_E * (2/r_p - 1/a_plus));

dv_p = v_p_plus - v_p_minus;
dv_flyby = norm(v_inf_plus - v_inf_minus);

%% === part 5 ===
% plot the two planetocentric hyperbolic arcs
figure("Name", "Planetocentric Powered Earth Fly-By"); hold on;

% --- set up ode solver conditions ---
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);
steps = 100;
tspan = linspace(0, 4800, steps);

% --- propagate and plot incoming planetocentric arc ---
y = [[r_p; 0; 0] [0; v_p_minus; 0;]];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), -tspan, y, options);
r_minus = Y(:, 1:3);
plot3(r_minus(:, 1), r_minus(:, 2), r_minus(:, 3), "r")

% --- propagate and plot outgoing planetocentric arc ---
y = [[r_p; 0; 0] [0; v_p_plus; 0;]];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y, options);
r_plus = Y(:, 1:3);
plot3(r_plus(:, 1), r_plus(:, 2), r_plus(:, 3), "g")

% --- plot earth surface texture ---
earth_img = imread("EarthTexture.jpg");
[x_earth, y_earth, z_earth] = sphere(50);
x_earth = r_E * x_earth; y_earth = r_E * y_earth; z_earth = -r_E * z_earth;
surface(x_earth, y_earth, z_earth, "FaceColor", "texturemap", ...
    "CData", earth_img, "EdgeColor", "none")

% --- finish up plot properties ---
legend("Incoming trajectory", "Outgoing trajectory", "Earth Texture")
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Powered Earth Fly-by Trajectory");
grid on; axis equal; hold off;

%% === part 5 (additional) ===
% i want to plot the heliocentric trajectory too
figure("Name", "Heliocentric Powered Earth Fly-By"); hold on;
steps = 1000;
tspan_helio = linspace(0, 60*60*24*30, steps);

% --- plot position of earth during fly-by ---
y = [R_E V_E];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), -tspan_helio, y, options);
R_E_minus = Y(:, 1:3);
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), 9*tspan_helio, y, options);
R_E_plus = Y(:, 1:3);
plot3(R_E_minus(:, 1), R_E_minus(:, 2), R_E_minus(:, 3), "b--");
plot3(R_E_plus(:, 1), R_E_plus(:, 2), R_E_plus(:, 3), "b");
scatter3(R_E(1), R_E(2), R_E(3), "filled", "g");

% --- propagate and plot incoming heliocentric arc ---
y = [R_E V_minus];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), -5*tspan_helio, y, options);
R_minus = Y(:, 1:3);
plot3(R_minus(:, 1), R_minus(:, 2), R_minus(:, 3), "r")

% --- propagate and plot outgoing heliocentric arc ---
y = [R_E V_plus];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), 5*tspan_helio, y, options);
R_plus = Y(:, 1:3);
plot3(R_plus(:, 1), R_plus(:, 2), R_plus(:, 3), "g")

% --- plot sun ---
scatter3(0, 0, 0, "filled", "y");

% --- finish up plot properties ---
legend("", "", "Earth", "Incoming Trajectory", "Outgoing Trajectory", "Sun")
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Powered Earth Fly-by Trajectory");
grid on; axis equal; hold off;

%% RESULTS DISPLAY
disp("     === LAB 4 EX. 2 RESULTS ===")
fprintf("       --- COMMON RESULTS ---\n")
fprintf("      v_inf_minus: [%.4f %.4f %.4f] km s^-1\n", v_inf_minus);
fprintf("       v_inf_plus: [%.4f %.4f %.4f] km s^-1\n", v_inf_plus);
fprintf("       turn angle: %.4f deg\n", rad2deg(turn_angle));
fprintf("pericentre radius: %.4f km\n", r_p);
fprintf("    fly-by height: %.4f km\n", h_ga);
fprintf("        delta-v p: %.4f km s^-1\n", dv_p);
fprintf("   fly-by delta-v: %.4f km s^-1\n", dv_flyby);
fprintf("       ecc before: %.4f\n", ecc_minus);
fprintf("       SMA before: %.4f km\n", a_minus);
fprintf("        ecc after: %.4f\n", ecc_plus);
fprintf("        SMA after: %.4f km\n", a_plus);
