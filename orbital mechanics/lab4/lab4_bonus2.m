%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
labs_d = cd(1:backs(end)); addpath([labs_d '\student_functions']); 
addpath([labs_d '\lib']); addpath([labs_d '\lib' '\timeConversion']);
clear; close all; clc;

%% Objective
% Achieve the same exit velocity as part 1 but by carrying out a powered
% trailing gravity assist

%% Data
mu_E    = astroConstants(13);   % earth gravitational parameter [km^2 s^-3]
mu_sun  = astroConstants(4);    % sun gravitational parameter [km^2 s^-3]
AU      = astroConstants(2);    % astronomical unit [km]

v_inf_plus = [10.0856; -6.4936; 0;];
v_inf_minus = [5.6934; 0; 0;];
R_E = [-1; 0; 0;] * AU;     % earth orbit radius [km]

%% Question 1
% Compute the heliocentric entry velocity [km s^-1]
V_E = sqrt(mu_sun / norm(R_E)) * [-R_E(2); R_E(1); 0;] / AU; % earth heliocentric velocity [km s^-1]
V_minus = V_E + v_inf_minus;

%% Question 2
% Compute the entry velocity in the SOI v_inf_minus [km s^-1]
u = cross(-impact_param, v_inf_plus); % vector normal to plane of hyperbola
u_hat = u / norm(u); % unit vector normal to plane of hyperbola
v_inf_minus = rodrigues(v_inf_plus, u_hat, turn_angle);

%% Question 3
% Compute the heliocentric velocity of the spacecraft before the flyby
% V_minus [km s^-1]
V_E = sqrt(mu_sun / norm(R_E)) * [-R_E(2); R_E(1); 0;] / AU; % earth heliocentric velocity [km s^-1]
V_minus = V_E + v_inf_minus;

%% RESULTS OUTPUT
disp("         === LAB 4 BONUS 2 RESULTS ===")
fprintf("1. impact parameter: [%.4f %.4f %.4f] km\n", impact_param);
fprintf("2.      v_inf_minus: [%.4f %.4f %.4f] km s^-1\n", v_inf_minus);
fprintf("3.          V_minus: [%.4f %.4f %.4f] km s^-1\n", V_minus);
