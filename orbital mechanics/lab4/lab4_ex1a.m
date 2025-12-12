%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
labs_d = cd(1:backs(end)); addpath([labs_d '\student_functions']); 
addpath([labs_d '\lib']); addpath([labs_d '\lib' '\timeConversion']);
clear; close all; clc;

% Exercise 1a
% Design a fly-by around the Earth for a fixed impact parameter and 
% different locations of the incoming asymptote.
%% 1. solve 2d hyperbola
% r_soi = r_p * (m_p / m_sun) ^ (2/5);

% --- constants ---
mu_E = astroConstants(13); % [km^3 s^-2]
mu_sun = astroConstants(4); % [km^3 s^-2]
AU = astroConstants(2); % [km]
R_E = AU .* [1; 0; 0;]; % [km]

% --- initial conditions ---
v_inf_minus = [15.1; 0; 0;]; % velocity before fly-by [km s^-1]
v_inf_mag = norm(v_inf_minus); % magnitude of velocity (const.) [km s^-1]
impact_param_mag = 9200; % impact parameter [km]
V_E = sqrt(mu_sun / R_E); % earth heliocentric velocity [km s^-1]

a = - mu_E / v_inf_mag^2;
turn_angle = 2 * atan(- a / impact_param_mag);

% build coordinate system with 
impact_param = [0; impact_param_mag; 0;];


%% 2. compute v_inf_plus for three asymptote positions
% 2.1 in front of the planet (decreased heliocentric velocity)
u = cross(-impact_param, v_inf_minus); % vector normal to plane of hyperbola
u_hat = u / norm(u); % unit vector normal to plane of hyperbola

v_inf_plus = rodrigues(v_inf_minus, u_hat, turn_angle);

% 2.2 
