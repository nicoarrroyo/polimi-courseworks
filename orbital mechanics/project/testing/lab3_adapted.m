%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
labs_d = cd(1:backs(end)); addpath([labs_d '\student_functions']); 
addpath([labs_d '\lib']); addpath([labs_d '\lib' '\timeConversion']);
%clear; close all; clc;

%% 1. Implement a function to compute the 𝚫𝒗𝐭𝐨𝐭(t1, t2)
% see transfer_cost function in student_functions directory

%% 2. Evaluate Δ𝑣tot for a grid of departure and arrival times covering ...
% the time windows provided for the Mercury-Earth leg.
% --- create array of departure and arrival times ---
dep_id = 1; arr_id = 3; steps = 100; mu_sun = astroConstants(4);
early_depart_mjd2000 = date2mjd2000([2040, 01, 01, 0, 00, 00]);
early_arrive_mjd2000 = date2mjd2000([2040, 04, 01, 0, 00, 00]);
late_depart_mjd2000 = date2mjd2000([2040, 03, 01, 0, 00, 00]);
late_arrive_mjd2000 = date2mjd2000([2040, 06, 01, 0, 00, 00]);

depart_times = linspace(early_depart_mjd2000, late_depart_mjd2000, steps);
arrive_times = linspace(early_arrive_mjd2000, late_arrive_mjd2000, steps);

% --- pre-calculate planet states ---
R_dep_list = zeros(steps, 3);
V_dep_list = zeros(steps, 3);
R_arr_list = zeros(steps, 3);
V_arr_list = zeros(steps, 3);

for i = 1:steps
    [R_dep_list(i, :), V_dep_list(i, :)] = ...
        get_planet_state(depart_times(i), dep_id, mu_sun);
    [R_arr_list(i, :), V_arr_list(i, :)] = ...
        get_planet_state(arrive_times(i), arr_id, mu_sun);
end

% --- conduct grid search ---
dvtot = NaN(steps, steps);
tof_days = NaN(steps, steps);
disp("conducting grid search"); tic
parfor i = 1:steps
    r1 = R_dep_list(i, :);
    v1 = V_dep_list(i, :);
    t1 = depart_times(i) * 24 * 3600;

    dv_row = NaN(1, steps);
    tof_row = NaN(1, steps);
    for j = 1:steps
        t2 = arrive_times(j) * 24 * 3600;
        tof = t2 - t1;
        
        if tof > 0
            r2 = R_arr_list(j, :);
            v2 = V_arr_list(j, :);

            [~, ~, ~, ERROR, v_t1, v_t2, ~, ~] = ...
                lambertMR(r1, r2, tof, mu_sun, 0, 0, 0, 0);
            if ERROR == 0
                dv_row(j) = norm(v_t1 - v1) + norm(v_t2 - v2);
                tof_row(j) = t2 - t1;
            end
        end
    end
    dvtot(i, :) = dv_row;
    tof_days(i, :) = tof_row / (24 * 3600);
end
disp("complete!"); toc

%% 3. Draw the porkchop plot of the Mercury-Earth Leg
figure("Name", "Mercury-Earth Leg Porkchop Plot"); hold on; grid on; axis equal;

% --- plot dv ---
dv_min = 10; dv_max = 30;
v_levels = dv_min : 1 : dv_max;
contour(depart_times, arrive_times, dvtot', v_levels, "ShowText", "off");
clim([dv_min, dv_max]);

% --- plot the constant time of flight lines ---
tof_levels = 50 : 50 : 200;
contour(depart_times, arrive_times, tof_days, tof_levels, "ShowText", "on");

% --- plot ---
colorbar;
xlabel("Departure Time (MJD2000)");
ylabel("Arrival Time (MJD2000)");
title("Porkchop Plot: Δv_{tot} for Mercury-Earth Leg");
hold off;

%% 4. Find the cheapest mission in terms of Δv
[min_val, min_idx_linear] = min(dvtot(:));
[row_idx, col_idx] = ind2sub(size(dvtot), min_idx_linear);

mjd_opt_dep_grid = depart_times(row_idx);
date_d = mjd20002date(mjd_opt_dep_grid);

mjd_opt_arr_grid = arrive_times(col_idx);
date_a = mjd20002date(mjd_opt_arr_grid);

tof_days_opt_grid = mjd_opt_arr_grid - mjd_opt_dep_grid;

% --- results output ---
fprintf("\n=== GRID SEARCH RESULTS ===\n");
fprintf("Min Delta V: %.4f km s^-1\n", min_val);
fprintf("Departure:   MJD2000 %.2f\n", mjd_opt_dep_grid);
fprintf("             Date    %.0f %.0f %.0f %.0f %.0f %.0f\n", date_d);
fprintf("Arrival:     MJD2000 %.2f\n", mjd_opt_arr_grid);
fprintf("             Date    %.0f %.0f %.0f %.0f %.0f %.0f\n", date_a);
fprintf("Total ToF (days):    %.2f\n", tof_days_opt_grid);

%% 5. Plot the transfer trajectory for this mission
% --- set up times and options ---
t1 = mjd_opt_dep_grid * 24 * 3600; % departure second
t2 = mjd_opt_arr_grid * 24 * 3600; % arrival second
tof = tof_days_opt_grid * 24 * 3600; % time of flight in seconds
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% --- get planet positions ---
[RD1, VD1] = get_planet_state(mjd_opt_dep_grid, dep_id, mu_sun);
[RA1, VA1] = get_planet_state(mjd_opt_dep_grid, arr_id, mu_sun);
[RD2, VD2] = get_planet_state(mjd_opt_arr_grid, dep_id, mu_sun);
[RA2, VA2] = get_planet_state(mjd_opt_arr_grid, arr_id, mu_sun);

[~, ~, ~, ~, V_transfer1, V_transfer2, ~, ~] = ...
    lambertMR( RD1, RA2, tof, mu_sun, 0, 0, 0, 0 );

% a. propagate the transfer arc
y_transfer = [RD1 V_transfer1]; % transfer arc initial state vector
tspan = linspace(t1, t2);

[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan, y_transfer, options);
R_transfer = Y(:, 1:3);
V_transfer = Y(:, 4:6);

% b. propagate the unused segment of the transfer arc
y_transfer_unused = [RA2 V_transfer2]; % transfer arc initial state vector
a_transfer = 1 / (2 / norm(RA2) - dot(V_transfer2, V_transfer2) / mu_sun); % semi major axis [km]
T_unused = 2*pi * sqrt(a_transfer^3 / mu_sun); % orbital period [s]
tspan_unused = linspace(t2, t2 + T_unused - t1); % integration time span array

[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_unused, y_transfer_unused, options);
R_transfer_unused = Y(:, 1:3);
V_transfer_unused = Y(:, 4:6);

% c. propagate the departure planet orbit
% i. during transfer
y_D = [RD1, VD1];
a_D = 1 / (2 / norm(RD1) - dot(VD1, VD1) / mu_sun); % semi major axis [km]
T_D = 2*pi * sqrt(a_D^3 / mu_sun); % orbital period [s]
tspan_D = linspace(t1, t2); % integration time span array

[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_D, y_D, options);
R_D = Y(:, 1:3);
V_D = Y(:, 4:6);
% ii. full orbit
% iii. transfer window

% d. propagate the arrival planet orbit
% i. during transfer
y_A = [RA1, VA1];
a_A = 1 / (2 / norm(RA1) - dot(VA1, VA1) / mu_sun); % semi major axis [km]
T_A = 2*pi * sqrt(a_A^3 / mu_sun); % orbital period [s]
tspan_A = linspace(t1, t2); % integration time span array

[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_A, y_A, options);
R_A = Y(:, 1:3);
V_A = Y(:, 4:6);
% ii. full orbit
% iii. transfer window

% --- plot ---
figure("Name", "Orbit Plot");

% a. transfer arc
plot3(R_transfer(:, 1), R_transfer(:, 2), R_transfer(:, 3), "y"); hold on;

% b. unused segment of transfer arc
%plot3(R_transfer_unused(:, 1), R_transfer_unused(:, 2), R_transfer_unused(:, 3), "--y");

% c. departure planet orbit
plot3(R_D(:, 1), R_D(:, 2), R_D(:, 3), "r"); % during transfer

% d. arrival planet orbit
plot3(R_A(:, 1), R_A(:, 2), R_A(:, 3), "g"); % during transfer

% e. initial departure planet position
scatter3(R_D(1, 1), R_D(1, 2), R_D(1, 3), "filled", "MarkerFaceColor", "r");
scatter3(R_D(end, 1), R_D(end, 2), R_D(end, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "r");

% f. final arrival planet position
scatter3(R_A(1, 1), R_A(1, 2), R_A(1, 3), "filled", "MarkerFaceColor", "g");
scatter3(R_A(end, 1), R_A(end, 2), R_A(end, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "g");

% g. sun position
scatter3(0, 0, 0, "filled", "MarkerFaceColor", "y");

% --- plot properties ---
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Two-body problem orbit");
legend(...
    "transfer arc", ...
    ..."unused transfer arc segment", ...
    "", ...
    "", ...
    "departure planet at departure", ...
    "departure planet at arrival", ...
    "arrival planet at arrival", ...
    "arrival planet at departure", ...
    "");
axis equal; grid on;
hold off;
