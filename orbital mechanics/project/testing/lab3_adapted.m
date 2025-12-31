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
leg1.dep_id = 1; leg1.arr_id = 3; steps = 100; mu_sun = astroConstants(4);
leg1.early_dept_mjd2000 = date2mjd2000([2040, 01, 01, 0, 00, 00]);
leg1.early_arr_mjd2000 = date2mjd2000([2040, 04, 01, 0, 00, 00]);
leg1.late_dept_mjd2000 = date2mjd2000([2040, 03, 01, 0, 00, 00]);
leg1.late_arr_mjd2000 = date2mjd2000([2040, 06, 01, 0, 00, 00]);

leg1.dep_times = linspace(leg1.early_dept_mjd2000, leg1.late_dept_mjd2000, steps);
leg1.arr_times = linspace(leg1.early_arr_mjd2000, leg1.late_arr_mjd2000, steps);

% --- pre-calculate planet states ---
leg1.R_dep_list = zeros(steps, 3);
leg1.V_dep_list = zeros(steps, 3);
leg1.R_arr_list = zeros(steps, 3);
leg1.V_arr_list = zeros(steps, 3);

for i = 1:steps
    [leg1.R_dep_list(i, :), leg1.V_dep_list(i, :)] = ...
        get_planet_state(leg1.dep_times(i), leg1.dep_id, mu_sun);
    [leg1.R_arr_list(i, :), leg1.V_arr_list(i, :)] = ...
        get_planet_state(leg1.arr_times(i), leg1.arr_id, mu_sun);
end

% --- conduct grid search ---
leg1.dvtot = NaN(steps, steps);
leg1.tof = NaN(steps, steps);
disp("conducting grid search"); tic
for i = 1:steps
    r1 = leg1.R_dep_list(i, :);
    v1 = leg1.V_dep_list(i, :);
    t1 = leg1.dep_times(i) * 24 * 3600;

    dv_row = NaN(1, steps);
    tof_row = NaN(1, steps);
    for j = 1:steps
        t2 = leg1.arr_times(j) * 24 * 3600;
        tof = t2 - t1;
        
        if tof > 0
            r2 = leg1.R_arr_list(j, :);
            v2 = leg1.V_arr_list(j, :);

            [~, ~, ~, ERROR, v_t1, v_t2, ~, ~] = ...
                lambertMR(r1, r2, tof, mu_sun, 0, 0, 0, 0);
            if ERROR == 0
                dv_row(j) = norm(v_t1 - v1) + norm(v_t2 - v2);
                tof_row(j) = t2 - t1;
            end
        end
    end
    leg1.dvtot(i, :) = dv_row;
    leg1.tof(i, :) = tof_row / (24 * 3600);
end
disp("complete!"); toc

%% 3. Draw the porkchop plot of the Mercury-Earth Leg
figure("Name", "Mercury-Earth Leg Porkchop Plot"); hold on; grid on; axis equal;

% --- plot dv ---
leg1.dv_min = 10; leg1.dv_max = 30;
leg1.v_levels = leg1.dv_min : 1 : leg1.dv_max;
contour(leg1.dep_times, leg1.arr_times, leg1.dvtot', leg1.v_levels, "ShowText", "off");
clim([leg1.dv_min, leg1.dv_max]);

% --- plot the constant time of flight lines ---
leg1.tof_levels = 50 : 50 : 200;
contour(leg1.dep_times, leg1.arr_times, leg1.tof, leg1.tof_levels, "ShowText", "on");

% --- plot ---
colorbar;
xlabel("Departure Time (MJD2000)");
ylabel("Arrival Time (MJD2000)");
title("Porkchop Plot: Δv_{tot} for Mercury-Earth Leg");
hold off;

%% 2. Evaluate Δ𝑣tot for a grid of departure and arrival times covering ...
% the time windows provided for the Earth-asteroid leg.
% --- create array of departure and arrival times ---
leg2.dep_id = leg1.arr_id; leg2.arr_id = 316801;
leg2.early_dept_mjd2000 = leg1.early_arr_mjd2000;
leg2.early_arr_mjd2000 = leg1.early_arr_mjd2000;
leg2.late_dept_mjd2000 = leg1.late_arr_mjd2000;
leg2.late_arr_mjd2000 = leg1.late_arr_mjd2000;

leg2.dep_times = linspace(leg2.early_dept_mjd2000, leg2.late_dept_mjd2000, steps);
leg2.arr_times = linspace(leg2.early_arr_mjd2000, leg2.late_arr_mjd2000, steps);

% --- pre-calculate planet states ---
leg2.R_dep_list = zeros(steps, 3);
leg2.V_dep_list = zeros(steps, 3);
leg2.R_arr_list = zeros(steps, 3);
leg2.V_arr_list = zeros(steps, 3);

for i = 1:steps
    [leg2.R_dep_list(i, :), leg2.V_dep_list(i, :)] = ...
        get_planet_state(leg2.dep_times(i), leg2.dep_id, mu_sun);
    [leg2.R_arr_list(i, :), leg2.V_arr_list(i, :)] = ...
        ephAsteroids(leg2.dep_times(i), leg2.arr_id);
end
[asteroid_arr.kep1, ~, asteroid_arr.M1] = ephAsteroids(leg2.dep_times(i), leg2.arr_id);

% --- conduct grid search ---
leg2.dvtot = NaN(steps, steps);
leg2.tof = NaN(steps, steps);
disp("conducting grid search"); tic
for i = 1:steps
    r1 = leg2.R_dep_list(i, :);
    v1 = leg2.V_dep_list(i, :);
    t1 = leg2.dep_times(i) * 24 * 3600;

    dv_row = NaN(1, steps);
    tof_row = NaN(1, steps);
    for j = 1:steps
        t2 = leg2.arr_times(j) * 24 * 3600;
        tof = t2 - t1;
        
        if tof > 0
            r2 = leg2.R_arr_list(j, :);
            v2 = leg2.V_arr_list(j, :);

            [~, ~, ~, ERROR, v_t1, v_t2, ~, ~] = ...
                lambertMR(r1, r2, tof, mu_sun, 0, 0, 0, 0);
            if ERROR == 0
                dv_row(j) = norm(v_t1 - v1) + norm(v_t2 - v2);
                tof_row(j) = t2 - t1;
            end
        end
    end
    leg2.dvtot(i, :) = dv_row;
    leg2.tof(i, :) = tof_row / (24 * 3600);
end
disp("complete!"); toc

%% 3. Draw the porkchop plot of the Mercury-Earth Leg
figure("Name", "Mercury-Earth Leg Porkchop Plot"); hold on; grid on; axis equal;

% --- plot dv ---
leg2.dv_min = 10; leg2.dv_max = 30;
leg2.v_levels = leg2.dv_min : 1 : leg2.dv_max;
contour(leg2.dep_times, leg2.arr_times, leg2.dvtot', leg2.v_levels, "ShowText", "off");
clim([leg2.dv_min, leg2.dv_max]);

% --- plot the constant time of flight lines ---
leg2.tof_levels = 50 : 50 : 200;
contour(leg2.dep_times, leg2.arr_times, leg2.tof, leg2.tof_levels, "ShowText", "on");

% --- plot ---
colorbar;
xlabel("Departure Time (MJD2000)");
ylabel("Arrival Time (MJD2000)");
title("Porkchop Plot: Δv_{tot} for Mercury-Earth Leg");
hold off;

%% 4. Find the cheapest mission in terms of Δv
[min_val, min_idx_linear] = min(leg1.dvtot(:));
[row_idx, col_idx] = ind2sub(size(leg1.dvtot), min_idx_linear);

mjd_opt_dep_grid = leg1.dep_times(row_idx);
date_d = mjd20002date(mjd_opt_dep_grid);

mjd_opt_arr_grid = leg1.arr_times(col_idx);
date_a = mjd20002date(mjd_opt_arr_grid);

leg1.tof_opt_grid = mjd_opt_arr_grid - mjd_opt_dep_grid;

% --- results output ---
fprintf("\n=== GRID SEARCH RESULTS ===\n");
fprintf("Min Delta V: %.4f km s^-1\n", min_val);
fprintf("Departure:   MJD2000 %.2f\n", mjd_opt_dep_grid);
fprintf("             Date    %.0f %.0f %.0f %.0f %.0f %.0f\n", date_d);
fprintf("Arrival:     MJD2000 %.2f\n", mjd_opt_arr_grid);
fprintf("             Date    %.0f %.0f %.0f %.0f %.0f %.0f\n", date_a);
fprintf("Total ToF (days):    %.2f\n", leg1.tof_opt_grid);

