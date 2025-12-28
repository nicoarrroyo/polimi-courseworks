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
a0 = 7200;          % Semi-major axis [km]
e0 = 0.0001;        % Eccentricity
i0 = deg2rad(95);   % Inclination [rad]
Omega0 = 0;         % RAAN [rad]
omega0 = 0;         % Argument of periapsis [rad]
TA0 = 0;            % True anomaly [rad]

kep0 = [a0; e0; i0; Omega0; omega0; TA0;];

% --- set intergration time: tspan ---
n_orbits = 500;
step_size = 5; % [s]

period = 2*pi * sqrt(a0^3 / params.mu);
tfinal = n_orbits * period;
n_steps = round(tfinal / step_size, 0);

tspan = linspace(t0, tfinal, n_steps);

% --- set ODE solver options: options ---
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

%% 2a. Propagate the given orbit using Gauss planetary equations
tic; disp("Gauss propogation...");
[T, Y_gauss] = ode113( ...
    @(t, k) eq_motion_gauss(t, k, @acc_pert_j2, params), ...
    tspan, ...
    kep0, ...
    options);
Y_gauss(:, 3:6) = rad2deg(Y_gauss(:, 3:6));
toc

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

tic; disp("Cartesian propogation");
[~, Y_car] = ode113( ...
    @(t,y) ode_2bp_j2(t, y, params.mu, params.J2, params.r_planet), ...
    tspan, ...
    y0, ...
    options);
r = Y_car(:, 1:3);
v = Y_car(:, 4:6);

[a, e, i, Omega, omega, TA] = car2kep(r, v, params.mu);
Y_car_kep_elements = [a, e, i, Omega, omega, TA];
Y_car_kep_elements(:, 3:6) = rad2deg(Y_car_kep_elements(:, 3:6));
Y_car_kep_elements(:, 4:6) = unwrap(Y_car_kep_elements(:, 4:6));
toc

% --- compute the error between the gaussian and cartesian solutions ---
error = zeros(size(Y_gauss));
if length(Y_gauss) == length(Y_car_kep_elements)
    error = (Y_gauss - Y_car_kep_elements) ./ Y_gauss;
    bad_sizes = 0;
else
    disp("arrays do not match in size");
    bad_sizes = 1;
end

%% 2c. Plot some characteristics for each element
% see later

%% 3a. Choose an appropriate cut-off period to remove oscillations
% The primary short-period oscillation due to J2 corresponds to the orbital 
% period[cite: 327]. We set the window width to one orbital period.
% The number of steps in one period is defined by the period divided by step_size.
window_width = round(period / step_size); 

%% 3b. Filter all elements
Y_gauss_unwrapped = Y_gauss;
% Convert angles (columns 3-6) to radians, unwrap, and convert back to degrees
Y_gauss_unwrapped(:, 3:6) = rad2deg(unwrap(deg2rad(Y_gauss(:, 3:6))));

% Apply the moving mean filter
Y_gauss_filtered = movmean(Y_gauss_unwrapped, window_width, 1);

%% 3c. Plot together the filtered and unfiltered results for each element
T_plot = T / period;

if ~bad_sizes
for kep_el = 1:width(Y_gauss)
    switch kep_el
        case 1
            fig_name = "Semi-major Axis";
            el_name = "a"; el_unit = "km";
        case 2
            fig_name = "Eccentricity";
            el_name = "e"; el_unit = "-";
        case 3
            fig_name = "Inclination";
            el_name = "i"; el_unit = "deg";
        case 4
            fig_name = "Right Ascension of Ascending Node";
            el_name = "Omega"; el_unit = "deg";
        case 5
            fig_name = "Argument of Perigee";
            el_name = "omega"; el_unit = "deg";
        case 6
            fig_name = "True Anomaly (TA)";
            el_name = "TA"; el_unit = "deg";
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
    plot(T_plot, Y_gauss_filtered(:, kep_el));
    title("All Orbits"); legend("Gaussian", "Cartesian", "Secular");
    xlabel("Time [n orbits]"); ylabel(el_name + " [" + el_unit + "]");
    xlim([T_plot(1), T_plot(end)]);
    hold off;

    subplot(3, 1, 3);
    portion_of_data = round(length(T_plot) / 10, 0);
    plot(T_plot(1:portion_of_data), Y_gauss(1:portion_of_data,kep_el)); hold on;
    plot(T_plot(1:portion_of_data), Y_car_kep_elements(1:portion_of_data,kep_el));
    plot(T_plot(1:portion_of_data), Y_gauss_filtered(1:portion_of_data, kep_el));
    title("10% of Orbits"); legend("Gaussian", "Cartesian", "Secular");
    xlabel("Time [n orbits]"); ylabel(el_name + " [" + el_unit + "]");
    hold off;
    
    % --- plot the filtered and unfiltered results ---
    % see subplots 3,1,2 and 3,1,3
end
end

%% 3d. Compare the slopes of the filtered Omega and omega with analytical J2
% Analytical approximations for secular rates 
% Constants
mu_E = params.mu;
R_E = params.r_planet;
J2 = params.J2;

% Common term in the J2 secular equations: (3/2)*n*J2*(Re/p)^2
common_factor = (1.5 * sqrt(mu_E) * J2 * R_E^2) / ((1 - e0^2)^2 * a0^3.5);

% Analytical Rates [rad/s]
dOmega_dt_anal = -common_factor * cos(i0); 
domega_dt_anal = common_factor * (2 - 2.5 * sin(i0)^2); % Equivalent to -[...](2.5sin^2 - 2)

% Convert analytical rates to [deg/s]
dOmega_dt_anal_deg = rad2deg(dOmega_dt_anal);
domega_dt_anal_deg = rad2deg(domega_dt_anal);

% Numerical Rates from Filtered Data [deg/s]
% use degree 1 polyfit on the filtered data to extract the slope.
coeffs_Omega = polyfit(T, Y_gauss_filtered(:, 4), 1);
dOmega_dt_num_deg = coeffs_Omega(1);

coeffs_omega = polyfit(T, Y_gauss_filtered(:, 5), 1);
domega_dt_num_deg = coeffs_omega(1);

% Display comparison
fprintf("--- Secular Evolution Comparison (J2) ---\n");
fprintf("RAAN (Omega) Rate [deg/s]:\n");
fprintf("  Analytical: %.5e\n", dOmega_dt_anal_deg);
fprintf("  Numerical:  %.5e\n", dOmega_dt_num_deg);
fprintf("  Rel Error:  %.2f %%\n", abs((dOmega_dt_anal_deg - dOmega_dt_num_deg)/dOmega_dt_anal_deg)*100);
fprintf("\n");
fprintf("Arg of Perigee (omega) Rate [deg/s]:\n");
fprintf("  Analytical: %.5e\n", domega_dt_anal_deg);
fprintf("  Numerical:  %.5e\n", domega_dt_num_deg);
fprintf("  Rel Error:  %.2f %%\n", abs((domega_dt_anal_deg - domega_dt_num_deg)/domega_dt_anal_deg)*100);
