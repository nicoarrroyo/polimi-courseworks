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
params.mu = astroConstants(13);
params.r_planet = astroConstants(23);
params.J2 = astroConstants(9);

% --- set initial time and state: t0, s0 ---
t0 = 0;
a0 = 7571;          % Semi-major axis [km]
e0 = 0.01;          % Eccentricity
i0 = deg2rad(87.9); % Inclination [deg]
Omega0 = pi;        % RAAN [deg]
omega0 = pi;        % Argument of periapsis [deg]
TA0 = 0;            % True anomaly [deg]

kep0 = [a0; e0; i0; Omega0; omega0; TA0;];

% --- set intergration time: tspan ---
n_orbits = 100;
step_size = 5; % [s]

period = 2*pi * sqrt(a0^3 / params.mu);
tfinal = n_orbits * period;
n_steps = round(tfinal / step_size, 0);

tspan = linspace(t0, tfinal, n_steps);

% --- set ODE solver options: options ---
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

%% 2a. Propagate the given orbit using Gauss planetary equations
[T, Y_gauss] = ode113( ...
    @(t, k) eq_motion_gauss(t, k, @acc_pert_j2, params), ...
    tspan, ...
    kep0, ...
    options);

%% 2b. Propagate the given orbit in Cartesian coordinates, and convert the 
% results to Keplerian elements
[r0, v0] = kep2car( ...          % state vector [km, km s^-1]
    a0, ...                      % a: semi major axis [km]
    e0, ...                      % e: eccentricity [-]
    i0, ...                      % i: inclination [rad]
    Omega0, ...                  % Omega: right ascension of the ascending node [rad]
    omega0, ...                  % omega: argument of perigee [rad]
    TA0, ...                     % TA: true anomaly [rad]
    params.mu, ...
    "radians" ...               % must specify unit of angles if not deg
    );
y0 = [r0; v0;];

[~, Y_car] = ode113( ...
    @(t,y) ode_2bp_j2(t, y, params.mu, params.J2, params.r_planet), ...
    tspan, ...
    y0, ...
    options);
r = Y_car(:, 1:3);
v = Y_car(:, 4:6);

[a, e, i, Omega, omega, TA] = car2kep(r, v, params.mu);
Y_car_kep_elements = [a, e, rad2deg(i), rad2deg(Omega), rad2deg(omega), rad2deg(TA)];

% --- compute the error between the gaussian and cartesian solutions ---
if length(Y_gauss) == length(Y_car_kep_elements)
    error = (Y_gauss - Y_car_kep_elements) ./ Y_gauss;
    bad_sizes = 0;
else
    disp("arrays do not match in size");
    bad_sizes = 1;
end

%% 2c. For each element
T_plot = T / period;

if ~bad_sizes
for kep_el = 1:width(Y_gauss)
    switch kep_el
        case 1
            fig_name = "Semi-major Axis";
            el_char = "a"; el_unit = "km";
        case 2
            fig_name = "Eccentricity";
            el_char = "e"; el_unit = "-";
        case 3
            fig_name = "Inclination";
            el_char = "i"; el_unit = "deg";
        case 4
            fig_name = "Right Ascension of Ascending Node";
            el_char = "Omega"; el_unit = "deg";
        case 5
            fig_name = "Argument of Perigee";
            el_char = "omega"; el_unit = "deg";
        case 6
            fig_name = "True Anomaly (TA)";
            el_char = "TA"; el_unit = "deg";
    end
    % --- plot the error between both solutions ---
    figure("Name", fig_name)
    subplot(3, 1, 1);
    plot(T_plot, error(:,kep_el));
    title("Cartesian vs Gaussian Error");
    xlabel("Time [n orbits]");
    ylabel("(" + el_name + "_g_a_u_s_s - " + el_name + "_c_a_r) / " + el_name + "_g_a_u_s_s");
    xlim([T_plot(1), T_plot(end)]);

    % --- plot both solutions together (gauss and cartesian) ---
    subplot(3, 1, 2);
    plot(T_plot, Y_gauss(:,kep_el)); hold on;
    plot(T_plot, Y_car_kep_elements(:,kep_el));
    title("All Orbits"); legend("Gaussian", "Cartesian");
    xlabel("Time [n orbits]"); ylabel(el_name + " [" + el_unit + "]");
    xlim([T_plot(1), T_plot(end)]);
    hold off;

    subplot(3, 1, 3);
    portion_of_data = round(length(T_plot) / 10, 0);
    plot(T_plot(1:portion_of_data), Y_gauss(1:portion_of_data,kep_el)); hold on;
    plot(T_plot(1:portion_of_data), Y_car_kep_elements(1:portion_of_data,kep_el));
    title("10% of Orbits"); legend("Gaussian", "Cartesian");
    xlabel("Time [n orbits]"); ylabel(el_name + " [" + el_unit + "]");
    hold off;
    
    % --- plot the filtered and unfiltered results ---
end
end
