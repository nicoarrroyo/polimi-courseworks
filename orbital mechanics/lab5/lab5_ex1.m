%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
labs_d = cd(1:backs(end)); addpath([labs_d '\student_functions']); 
addpath([labs_d '\lib']); addpath([labs_d '\lib' '\timeConversion']);
clear; close all; clc;

% Exercise: propagate an earth orbit perturbed by J2 using Gauss planetary
% equations and study the results
% 1. Implement the code for orbit propagation with Gauss planetary equations
%% --- a. implement a function with the equations of motion ---
% see eq_motion_gauss

%% --- b. implement a function for the perturbing acceleration due to J2
% see acc_pert_j2

%% --- c. implement a function for orbit propagation ---
% --- set physical parameters ---
params.mu_E = astroConstants(13);
params.r_E = astroConstants(23);
params.J2 = astroConstants(9);

% --- set initial time and state: t0, s0 ---
t0 = 0;
kep.a0 = 7571;          % Semi-major axis [km]
kep.e0 = 0.01;          % Eccentricity
kep.i0 = deg2rad(87.9); % Inclination [deg]
kep.Omega0 = pi;        % RAAN [deg]
kep.omega0 = pi;        % Argument of periapsis [deg]
kep.TA0 = 0;            % True anomaly [deg]

kep.s0 = [kep.a0; kep.e0; kep.i0; kep.Omega0; kep.omega0; kep.TA0;];

[car.r0, car.v0] = kep2car( ... % initial state vector [km, km s^-1]
    kep.a0, ...                 % a: semi major axis [km]
    kep.e0, ...                 % e: eccentricity [-]
    rad2deg(kep.i0), ...        % i: inclination [deg]
    rad2deg(kep.Omega0), ...    % Omega: right ascension of the ascending node [deg]
    rad2deg(kep.omega0), ...    % omega: argument of perigee [deg]
    rad2deg(kep.TA0));          % TA: true anomaly [deg]
car.r0 = car.r0'; car.v0 = car.v0'; % transpose to match Lambert function dimensions

car.s0 = [car.r0; car.v0;];

% --- set intergration time: tspan ---
n_orbits = 100;
step_size = 5; % [s]

T = 2*pi * sqrt(kep.a0^3 / mu_E);
tfinal = n_orbits * T;
n_steps = round(tfinal / step_size, 0);

tspan = linspace(t0, tfinal, n_steps);

% --- set ODE solver options: options ---
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% --- numerical integration of equations of motion ---
[T, S] = ode113( ...
    @(t, s) eq_motion( ...
    t, s, @(t, s) acc_pert_fun(t, s, parameters), parameters), ...
    tspan, s0, options);

% --- analyse and plot results ---














%% %% 1. Create a "wrapper" for the perturbations that already includes your constants
% % This makes acc_pert_fun look like it only needs (t, s)
% acc_wrap = @(t, s) acc_pert_fun(t, s, parameters);
% 
% % 2. Create a "wrapper" for the main equations of motion
% % We "bake in" the perturbation function and the physical parameters
% ode_fun = @(t, s) eq_motion(t, s, acc_wrap, parameters);
% 
% % 3. Pass the simple (t, s) function to the solver
% % Now the solver only sees two inputs, which is what it expects
% [T, S] = ode113(ode_fun, tspan, s0, options);

% %% 1. Wrap the perturbation model with its specific parameters
% pert_handle = @(t, s) acc_pert_fun(t, s, parameters);
% 
% % 2. Wrap the equations of motion with the perturbation handle and parameters
% ode_fun = @(t, s) eq_motion(t, s, pert_handle, parameters);
% 
% % 3. Call the solver using the clean "ode_fun" handle
% [T, S] = ode113(ode_fun, tspan, s0, options);
% 
% %% No @ symbols used here
% [T, S] = ode113("bridge_function", tspan, s0, options);
% 
% % No @ symbols used here
% [T, S] = ode113("bridge_function", tspan, s0, options);
% 
% function ds = bridge_function(t, s)
%     % You have to hard-code or load the parameters because 
%     % the solver doesn't pass them to this function!
%     load("constants.mat"); 
% 
%     % We have to manually define which perturbation to use
%     pert_vec = acc_pert_fun(t, s, parameters);
% 
%     % Now call the actual equations of motion
%     ds = eq_motion(t, s, pert_vec, parameters);
% end
