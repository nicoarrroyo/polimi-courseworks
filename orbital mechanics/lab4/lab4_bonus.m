%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
labs_d = cd(1:backs(end)); addpath([labs_d '\student_functions']); 
addpath([labs_d '\lib']); addpath([labs_d '\lib' '\timeConversion']);
clear; close all; clc;

%% === Objective 1 ===
% Design a leading flyby around the Earth, given the following data on the 
% turning angle and on the exit velocity in the SOI

%% === Data 1 ===
mu_E    = astroConstants(13);   % earth gravitational parameter [km^2 s^-3]
mu_sun  = astroConstants(4);    % sun gravitational parameter [km^2 s^-3]
AU      = astroConstants(2);    % astronomical unit [km]

turn_angle = deg2rad(25.8230);
v_inf_plus = [10.0856; -6.4936; 0;];
R_E = [0; -1; 0;] * AU;     % earth orbit radius [km]

%% === Question 1 ===
% Compute the value of the impact parameter [km]
a = -mu_E / (norm(v_inf_plus)^2);
impact_param_mag = -a / tan(turn_angle / 2);
impact_param = [0; impact_param_mag; 0;];

%% === Question 2 ===
% Compute the entry velocity in the SOI v_inf_minus [km s^-1]
u = cross(-impact_param, v_inf_plus); % vector normal to plane of hyperbola
u_hat = u / norm(u); % unit vector normal to plane of hyperbola
v_inf_minus = rodrigues(v_inf_plus, u_hat, turn_angle);

%% === Question 3 ===
% Compute the heliocentric velocity of the spacecraft before the flyby
% V_minus [km s^-1]
% --- earth heliocentric velocity [km s^-1] ---
V_E = sqrt(mu_sun / norm(R_E)) * [-R_E(2); R_E(1); 0;] / AU;

V_minus = V_E + v_inf_minus;

%% === RESULTS OUTPUT ===
fprintf("      === LAB 4 BONUS 1 RESULTS ===\n")
fprintf("1. impact parameter: [%.4f %.4f %.4f] km\n", impact_param);
fprintf("2.      v_inf_minus: [%.4f %.4f %.4f] km s^-1\n", v_inf_minus);
fprintf("3.          V_minus: [%.4f %.4f %.4f] km s^-1\n", V_minus);

%% === Objective 2 ===
% Achieve the same exit velocity as part 1 but by carrying out a powered
% trailing gravity assist

%% === Data 2 ===
r_E         = astroConstants(23);
v_inf_plus  = [10.0856; -6.4936; 0;];
v_inf_minus = [5.6934; 0; 0;];
R_E         = [-1; 0; 0;] * AU; % earth orbit radius [km]

%% === Question 4 ===
% Compute the heliocentric entry velocity [km s^-1]
% --- earth heliocentric velocity [km s^-1] ---
V_E = sqrt(mu_sun / norm(R_E)) * [-R_E(2); R_E(1); 0;] / AU;

V_minus = V_E + v_inf_minus;

%% === Question 5 ===
% Compute the turning angle of the new flyby [deg]
turn_angle = acos(dot(v_inf_minus, v_inf_plus) / (norm(v_inf_minus) * norm(v_inf_plus)));
turn_angle_deg = rad2deg(turn_angle);

%% === Question 6 ===
% Compute the semi-major axis of the outgoing arc a_plus [km]
a_plus = -mu_E / (norm(v_inf_plus)^2);

%% === Question 7 ===
% Compute the velocity impulse provided at pericentre dv_p [km_s]
% --- solve the non-linear system for r_p ---
eq = @(r_p) turn_angle - ...
    asin(1 / (1 + (r_p * norm(v_inf_plus)^2) / mu_E)) - ...
    asin(1 / (1 + (r_p * norm(v_inf_minus)^2) / mu_E));
r_p = fzero(eq, r_E + 500);

% --- check its validity ---
h_atm       = 500;          % height of earth atmosphere from sea-level [km]
r_p_crit    = r_E + h_atm;  % critical fly-by pericentre radius

if r_p < r_p_crit
    fprintf("\n!!! DANGER DANGER DANGER !!!\n")
    fprintf("\nfly-by pericentre radius might be too low!\n")
    fprintf("\nTERRAIN. PULL UP. TERRAIN. PULL UP.\n")
    fprintf("\n!!! DANGER DANGER DANGER !!!\n")
end

% --- compute pericentre velocity before and after manouvre ---
a_minus = -mu_E / (norm(v_inf_minus)^2); % 
v_p_minus = sqrt(mu_E * (2/r_p - 1/a_minus));
v_p_plus = sqrt(mu_E * (2/r_p - 1/a_plus));

dv_p = v_p_plus - v_p_minus;

%% === RESULTS OUTPUT ===
fprintf("\n      === LAB 4 BONUS 2 RESULTS ===\n")
fprintf("4.    V_minus: [%.4f %.4f %.4f] km s^-1\n", V_minus);
fprintf("5. turn angle: %.4f deg\n", turn_angle_deg);
fprintf("6.     a_plus: %.4f km\n", a_plus);
fprintf("7.       dv_p: [%.4f %.4f %.4f] km s^-1\n", dv_p);
