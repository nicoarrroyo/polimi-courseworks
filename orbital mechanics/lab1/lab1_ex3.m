%% configure paths
script_path = fileparts(mfilename("fullpath"));
backs = strfind(script_path, "\"); labs_dir = script_path(1:backs(end));
addpath([labs_dir '\student_functions']); addpath([labs_dir '\lib']);
clear; close all; clc;

%% Exercise 3 - Kepler's Equation
J2 = astroConstants(9); % second zonal harmonic
mu_E = astroConstants(13); % [km^3 s^-2]
R_e = astroConstants(23); % earth radius [km]

% initial conditions
E0 = 0; % initial eccentric anomaly [rad]
t0 = 0; % initial time [s]

% others
a = 7000; % semi major axis [km]
k = 2; % periods of the orbit
N = 200; % steps of orbit
e_values = [0 0.2 0.4 0.6 0.8 0.95]; % array of eccentricities
n = sqrt(mu_E / a^3); % mean motion [km s^-1]
period = 2 * pi / n; % time period [s]

tspan = linspace(0, k*period, N); % time span

% get results
E_results = zeros(length(e_values), N);
for i = 1:length(e_values)
    e = e_values(i);
    for j = 1:N
        E_results(i, j) = kepler_solver(tspan(j), e, a, mu_E, t0, E0);
    end
end

% plot results
% 2d plot
figure("Name", "Eccentric Anomaly")
hold on
for i = 1:height(E_results)
    plot(tspan/period, rad2deg(E_results(i, :)))
end
xlabel("Time [T]"); ylabel("E [deg]")
title("Eccentric Anomaly against Time for Different Eccentricities")
grid on;
legend("e=0", "e=0.2", "e=0.4", "e=0.6", "e=0.8", "e=0.95");
hold off

% surface plot
figure("Name", "Eccentric Anomaly vs Time vs Eccentricity")
[T_grid, E_grid] = meshgrid(tspan / period, e_values);
surf(T_grid, E_grid, rad2deg(E_results));
xlabel("Time [T]");
ylabel("Eccentricity e");
zlabel("Eccentric Anomaly E");
title("Eccentric Anomaly as a function of Time and Eccentricity");
colorbar;
