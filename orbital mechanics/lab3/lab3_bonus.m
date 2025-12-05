%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
labs_d = cd(1:backs(end)); addpath([labs_d '\student_functions']); 
addpath([labs_d '\lib']); addpath([labs_d '\lib' '\timeConversion']);
clear; close all; clc;

%% Question 1
% a. positions
r1 = [-1964.809 2821.834 596.808];
v1 = [-2.902 -2.044 0.107];

r2 = [4836.089 -6945.559 -1468.959];
v2 = [1.629 1.147 -0.061];

% b. lambert
tof = 1*3600 + 58.76*60;
mu_mars = astroConstants(14); % mars gravitational parameter [km^3/s^2]

[~, ~, ~, ~, v_transfer1, v_transfer2, ~, ~] = ...
    lambertMR( r1, r2, tof, mu_mars, 0, 0, 0, 0 );

% 3. Compute the total cost of the manoeuvre (Δv1 + Δv2)
dv1 = norm(abs(v_transfer1 - v1));
dv2 = norm(abs(v_transfer2 - v2));
dvtot = dv1 + dv2;

% 4. Results
fprintf("===== QUESTION 1 RESULTS =====\n");
fprintf("total Δv:         %.4f km s^-1\n", dvtot);

%% Question 2
[a, ~, ~, ~, ~, ~] = car2kep(r1, v_transfer1);

fprintf("\n===== QUESTION 2 RESULTS =====\n");
fprintf("semi-major axis:  %.4f km s^-1\n", a);

%% Question 3
% hohmann transfer

%% Question 4


