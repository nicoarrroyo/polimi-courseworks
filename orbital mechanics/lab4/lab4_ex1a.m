%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
labs_d = cd(1:backs(end)); addpath([labs_d '\student_functions']); 
addpath([labs_d '\lib']); addpath([labs_d '\lib' '\timeConversion']);
clear; close all; clc;

% Exercise 1a
% Design a fly-by around the Earth for a fixed impact parameter and 
% different locations of the incoming asymptote.
%% 1. solve 2d hyperbola
% --- constants ---
mu_E                = astroConstants(13);   % earth grav. param. [km^3 s^-2]
mu_sun              = astroConstants(4);    % sun grav. param. [km^3 s^-2]
AU                  = astroConstants(2);    % astronomical unit [km]
R_E                 = AU .* [0; 1; 0;];     % earth orbital radius [km]
r_E                 = astroConstants(23);   % earth mean radius [km]
impact_param_mag    = 9200;                 % impact parameter [km]

% --- initial conditions ---
v_inf_minus = [15.1; 0; 0;]; % velocity before fly-by [km s^-1]

% --- calculated terms ---
a           = - mu_E / (norm(v_inf_minus) ^ 2); % semi-major axis [km]
turn_angle  = 2 * atan(- a / impact_param_mag); % turn angle, delta [rad]
ecc         = 1 / sin(turn_angle / 2);          % eccentricity [-]
r_p         = a * (1 - ecc);                    % pericentre radius [km]

%% 2. compute v_inf_plus for three asymptote positions
% 2.1 in front of the planet (decreased heliocentric velocity)
impact_param_leading = [0; impact_param_mag; 0;];
u_leading = cross(-impact_param_leading, v_inf_minus); % vector normal to plane of hyperbola
u_hat_leading = u_leading / norm(u_leading); % unit vector normal to plane of hyperbola
v_inf_plus_leading = rodrigues(v_inf_minus, u_hat_leading, turn_angle);

% 2.2 behind the planet (increased heliocentric velocity)
impact_param_trailing = [0; impact_param_mag; 0;];
u_trailing = cross(impact_param_trailing, v_inf_minus);
u_hat_trailing = u_trailing / norm(u_trailing);
v_inf_plus_trailing = rodrigues(v_inf_minus, u_hat_trailing, turn_angle);

% 2.3 under the planet (change of plane)
impact_param_under = [0; 0; -impact_param_mag;];
u_under = cross(impact_param_under, v_inf_minus);
u_hat_under = u_under / norm(u_under);
v_inf_plus_under = rodrigues(v_inf_minus, u_hat_under, turn_angle);

%% 3. compute V_minus, V_plus, and the incoming and outgoing heliocentric arcs
% --- common ---
V_E = sqrt(mu_sun / R_E)'; % earth heliocentric velocity [km s^-1]
dv = abs(norm(v_inf_minus) - norm(v_inf_plus_under));

% --- set ode solver conditions ---
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);
tspan = linspace(0, 10000, 500);

% 3.1 leading-side fly-by
V_minus_leading = V_E + v_inf_minus;
V_plus_leading = V_E + v_inf_plus_leading;

y_leading = [[-5*r_E; impact_param_mag; 0;] v_inf_minus];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y_leading, options);
r_leading = Y(:, 1:3);

% 3.2 trailing-side fly-by
V_minus_trailing = V_E + v_inf_minus;
V_plus_trailing = V_E + v_inf_plus_trailing;

y_trailing = [[-5*r_E; -impact_param_mag; 0;] v_inf_minus];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y_trailing, options);
r_trailing = Y(:, 1:3);

% 3.3 under-the-planet fly-by
V_minus_under = V_E + v_inf_minus;
V_plus_under = V_E + v_inf_plus_under;

y_under = [[-5*r_E; 0; -impact_param_mag;] v_inf_minus];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y_under, options);
r_under = Y(:, 1:3);

%% results output
% --- create figure ---
figure("Name", "Leading, Trailing, Under"); hold on;

% --- plot leading, trailing, and under side fly-by trajectories
plot3(r_leading(:, 1), r_leading(:, 2), r_leading(:, 3), "r");
plot3(r_trailing(:, 1), r_trailing(:, 2), r_trailing(:, 3), "g");
plot3(r_under(:, 1), r_under(:, 2), r_under(:, 3), "b");

% --- plot earth surface texture ---
earth_img = imread("EarthTexture.jpg");
[x_earth, y_earth, z_earth] = sphere(50);
x_earth = r_E * x_earth; y_earth = r_E * y_earth; z_earth = -r_E * z_earth;
surface(x_earth, y_earth, z_earth, "FaceColor", "texturemap", ...
    "CData", earth_img, "EdgeColor", "none")

% --- finish up plot properties ---
legend("Leading-side", "Trailing-side", "Under-side")

xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Earth Fly-by Trajectories");

axis equal; grid on;
xlim([-2*r_E 16*r_E]);
hold off;

% --- output final results ---
disp("=== LAB 3 EX. 1a RESULTS ===")
fprintf("\n--- COMMON RESULTS ---\n")
fprintf("       turn angle: %.4f deg\n", rad2deg(turn_angle));
fprintf("pericentre radius: %.4f km s^-1\n", r_p);
fprintf("   fly-by delta-v: %.4f km s^-1\n", dv);
fprintf("  semi-major axis: %.4f km\n", a);
fprintf("     eccentricity: %.4f km\n", ecc);

fprintf("\n--- LEADING-SIDE FLY-BY ---\n")
fprintf("V\x207B:   %.4f %.4f %.4f km\n", V_minus_trailing);
fprintf("V\x207A:   %.4f %.4f %.4f km\n", V_plus_trailing);
fprintf("v\x207A\x221e:  %.4f %.4f %.4f km\n", v_inf_plus_leading);

fprintf("\n--- TRAILING-SIDE FLY-BY ---\n")
fprintf("V\x207B:   %.4f %.4f %.4f km\n", V_minus_trailing);
fprintf("V\x207A:   %.4f %.4f %.4f km\n", V_plus_trailing);
fprintf("v\x207A\x221e:  %.4f %.4f %.4f km\n", v_inf_plus_trailing);

fprintf("\n--- UNDER-SIDE FLY-BY ---\n")
fprintf("V\x207B:   %.4f %.4f %.4f km\n", V_minus_under);
fprintf("V\x207A:   %.4f %.4f %.4f km\n", V_plus_under);
fprintf("v\x207A\x221e:  %.4f %.4f %.4f km\n", v_inf_plus_under);
