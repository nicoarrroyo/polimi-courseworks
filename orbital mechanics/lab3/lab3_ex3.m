%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
labs_d = cd(1:backs(end)); addpath([labs_d '\student_functions']); 
addpath([labs_d '\lib']); addpath([labs_d '\lib' '\timeConversion']);
clear; close all; clc;

% Exercise 3: Mars Express
% Design an interplanetary transfer with minimum 𝚫𝒗𝐭𝐨𝐭 between Earth and 
% Mars, under the following mission requirements:
% time_depart_early = 2003 April 1; time_depart_late = 2003 August 1;
% time_arrive_early = 2003 September 1; time_arrive_late = 2004 March 1;
% No need to consider planetary insertion date.

%% 1. Implement a function to compute the 𝚫𝒗𝐭𝐨𝐭(t1, t2)
% see transfer_cost function in student_functions directory

%% 2. Evaluate Δ𝑣tot for a grid of departure and arrival times covering ...
% the time windows provided.
early_depart_mjd2000 = date2mjd2000([2003, 04, 01, 0, 00, 00]);
early_arrive_mjd2000 = date2mjd2000([2003, 09, 01, 0, 00, 00]);

late_depart_mjd2000 = date2mjd2000([2003, 08, 01, 0, 00, 00]);
late_arrive_mjd2000 = date2mjd2000([2004, 03, 01, 0, 00, 00]);

steps = 100;
depart_times = linspace(early_depart_mjd2000, late_depart_mjd2000, steps);
arrive_times = linspace(early_arrive_mjd2000, late_arrive_mjd2000, steps);

dvtot = zeros(length(depart_times), length(arrive_times));
tof_days = zeros(length(depart_times), length(arrive_times));

for i = 1:length(depart_times)
    for j = 1:length(arrive_times)
        dvtot(i, j) = transfer_cost(...
            depart_times(i), arrive_times(j), 3, 4);
        
        tof_days(i, j) = arrive_times(j) - depart_times(i);
    end
end

%% 3. Draw the porkchop plot of the Mars Express Mission
figure("Name", "Mars Express Porkchop Plot"); hold on; grid on; axis equal;

% plot dv
dv_min = 5; dv_max = 10;
v_levels = dv_min : 1 : dv_max;
contour(depart_times, arrive_times, dvtot', v_levels, "ShowText", "on");
clim([dv_min, dv_max]);

% plot the constant time of flight lines
tof_levels = 100 : 100 : 400;
contour(depart_times, arrive_times, tof_days, tof_levels, "ShowText", "on");

colorbar;
xlabel("Departure Time (MJD2000)");
ylabel("Arrival Time (MJD2000)");
title("Porkchop Plot: Δv_{tot} for Mars Express Mission");
hold off;

%% 4. Find the cheapest mission in terms of Δv
[min_val, min_idx_linear] = min(dvtot(:));
[row_idx, col_idx] = ind2sub(size(dvtot), min_idx_linear);

mjd_opt_dep_grid = depart_times(row_idx);
date_d = mjd20002date(mjd_opt_dep_grid);

mjd_opt_arr_grid = arrive_times(col_idx);
date_a = mjd20002date(mjd_opt_arr_grid);

tof_days_opt_grid = mjd_opt_arr_grid - mjd_opt_dep_grid;

% results output
fprintf("\n=== GRID SEARCH RESULTS ===\n");
fprintf("Min Delta V: %.4f km s^-1\n", min_val);
fprintf("Departure:   MJD2000 %.2f\n", mjd_opt_dep_grid);
fprintf("             Date    %.0f %.0f %.0f %.0f %.0f %.0f\n", date_d);
fprintf("Arrival:     MJD2000 %.2f\n", mjd_opt_arr_grid);
fprintf("             Date    %.0f %.0f %.0f %.0f %.0f %.0f\n", date_a);
fprintf("Total ToF (days):    %.2f\n", tof_days_opt_grid);

%% 5. Plot the transfer trajectory for this mission
% constants
mu_sun = astroConstants(4);
earth_id = 3; mars_id = 4;

% get initial planet positions
[RE1, VE1] = get_planet_state(mjd_opt_dep_grid, earth_id, mu_sun);
[RM1, VM1] = get_planet_state(mjd_opt_dep_grid, mars_id, mu_sun);

% get final planet positions
[RE2, VE2] = get_planet_state(mjd_opt_arr_grid, earth_id, mu_sun);
[RM2, VM2] = get_planet_state(mjd_opt_arr_grid, mars_id, mu_sun);

% convert days to seconds
seconds_per_day = 24 * 3600;
t1 = mjd_opt_dep_grid * seconds_per_day; % departure second
t2 = mjd_opt_arr_grid * seconds_per_day; % arrival second
tof = tof_days_opt_grid * seconds_per_day; % time of flight in seconds

[~, ~, ~, ~, V_transfer1, V_transfer2, ~, ~] = ...
    lambertMR( RE1, RM2, tof, mu_sun, 0, 0, 0, 0 );

% set ode solver conditions
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% a. propagate the transfer arc
y_transfer = [RE1 V_transfer1]; % transfer arc initial state vector
tspan = linspace(t1, t2);

[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan, y_transfer, options);
R_transfer = Y(:, 1:3);
V_transfer = Y(:, 4:6);

% b. propagate the unused segment of the transfer arc
y_transfer_unused = [RM2 V_transfer2]; % transfer arc initial state vector
a_transfer = 1 / (2 / norm(RM2) - dot(V_transfer2, V_transfer2) / mu_sun); % semi major axis [km]
T_unused = 2*pi * sqrt(a_transfer^3 / mu_sun); % orbital period [s]
tspan_unused = linspace(t2, t2 + T_unused - t1); % integration time span array

[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_unused, y_transfer_unused, options);
R_transfer_unused = Y(:, 1:3);
V_transfer_unused = Y(:, 4:6);

% c. propagate the earth orbit
y_earth = [RE1, VE1];
a_earth = 1 / (2 / norm(RE1) - dot(VE1, VE1) / mu_sun); % semi major axis [km]
T_earth = 2*pi * sqrt(a_earth^3 / mu_sun); % orbital period [s]
tspan_earth = linspace(t1, t2); % integration time span array

[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_earth, y_earth, options);
R_earth = Y(:, 1:3);
V_earth = Y(:, 4:6);
% ii. full orbit
% iii. transfer window

% d. propagate the mars orbit
% i. during transfer
y_mars = [RM1, VM1];
a_mars = 1 / (2 / norm(RM1) - dot(VM1, VM1) / mu_sun); % semi major axis [km]
T_mars = 2*pi * sqrt(a_mars^3 / mu_sun); % orbital period [s]
tspan_mars = linspace(t1, t2); % integration time span array

[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_mars, y_mars, options);
R_mars = Y(:, 1:3);
V_mars = Y(:, 4:6);
% ii. full orbit
% iii. transfer window

% plot
figure("Name", "Orbit Plot");

% a. transfer arc
plot3(R_transfer(:, 1), R_transfer(:, 2), R_transfer(:, 3), "y"); hold on;

% b. unused segment of transfer arc
plot3(R_transfer_unused(:, 1), R_transfer_unused(:, 2), R_transfer_unused(:, 3), "--y");

% c. earth orbit
plot3(R_earth(:, 1), R_earth(:, 2), R_earth(:, 3), "r"); % during transfer
%plot3(R_earth(:, 1), R_earth(:, 2), R_earth(:, 3), "r"); % full orbit

% d. mars orbit
plot3(R_mars(:, 1), R_mars(:, 2), R_mars(:, 3), "g"); % during transfer
%plot3(R_mars(:, 1), R_mars(:, 2), R_mars(:, 3), "g"); % full orbit

% e. initial earth position
scatter3(R_earth(1, 1), R_earth(1, 2), R_earth(1, 3), "filled", "MarkerFaceColor", "g");
scatter3(R_earth(end, 1), R_earth(end, 2), R_earth(end, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "g");

% f. final mars position
scatter3(R_mars(1, 1), R_mars(1, 2), R_mars(1, 3), "filled", "MarkerFaceColor", "r");
scatter3(R_mars(end, 1), R_mars(end, 2), R_mars(end, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "r");

% g. sun position
scatter3(0, 0, 0, "filled", "MarkerFaceColor", "y");

% plot properties
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Two-body problem orbit");
legend(...
    "transfer arc", ...
    "unused transfer arc segment", ...
    "", ...
    "", ...
    "earth at departure", ...
    "earth at arrival", ...
    "mars at arrival", ...
    "mars at departure", ...
    "");
axis equal; grid on;
hold off;

