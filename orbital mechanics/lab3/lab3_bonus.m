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
[a, ~, ~, ~, ~, ~] = car2kep(r1, v_transfer1, mu_mars);

fprintf("\n===== QUESTION 2 RESULTS =====\n");
fprintf("semi-major axis:  %.4f km s^-1\n", a);

%% Question 3
% hohmann transfer

%% Question 4
% 1. Evaluate Δ𝑣tot for a grid of departure and arrival times within the 
% given time windows.
departure_id = 1; arrival_id = 2;
early_depart_mjd2000 = date2mjd2000([2024, 06, 01, 0, 00, 00]);
early_arrive_mjd2000 = date2mjd2000([2024, 12, 01, 0, 00, 00]);

late_depart_mjd2000 = date2mjd2000([2026, 11, 01, 0, 00, 00]);
late_arrive_mjd2000 = date2mjd2000([2027, 06, 01, 0, 00, 00]);

steps = 5000;
depart_times = linspace(early_depart_mjd2000, late_depart_mjd2000, steps);
arrive_times = linspace(early_arrive_mjd2000, late_arrive_mjd2000, steps);

%% 4. Find the cheapest mission in terms of Δv (with launcher constraint)
v_inf = 7; % specific excess velocity [km s^-1]

dv1 = NaN(length(depart_times), length(arrive_times));
dv2 = NaN(length(depart_times), length(arrive_times));

for i = 1:length(depart_times)
    for j = 1:length(arrive_times)
        if arrive_times(j) < depart_times(i)
            continue
        end
        [dv1(i, j), dv2(i, j)] = constrained_transfer_cost(...
            depart_times(i), ...
            arrive_times(j), ...
            departure_id, ...
            arrival_id, ...
            v_inf);
    end
end

dvtot_constrained = dv1 + dv2;
[min_val_constrained, min_idx_linear_constrained] = min(dvtot_constrained(:));
[row_idx_constrained, col_idx_constrained] = ind2sub(size(dvtot_constrained), min_idx_linear_constrained);

mjd_opt_dep_grid_constrained = depart_times(row_idx_constrained);
date_d_constrained = mjd20002date(mjd_opt_dep_grid_constrained);

mjd_opt_arr_grid_constrained = arrive_times(col_idx_constrained);
date_a_constrained = mjd20002date(mjd_opt_arr_grid_constrained);

tof_days_opt_grid_constrained = mjd_opt_arr_grid_constrained - mjd_opt_dep_grid_constrained;

% results output
fprintf("\n=== QUESTION 4 RESULTS ===\n");
fprintf("Min Delta V:      %.4f km s^-1\n", min_val_constrained);
fprintf("Departure:        MJD2000 %.2f\n", mjd_opt_dep_grid_constrained);
fprintf("                  Date    %.0f %.0f %.0f %.0f %.0f %.2f\n", date_d_constrained);
fprintf("Arrival:          MJD2000 %.2f\n", mjd_opt_arr_grid_constrained);
fprintf("                  Date    %.0f %.0f %.0f %.0f %.0f %.2f\n", date_a_constrained);
fprintf("ToF (days):       %.2f\n", tof_days_opt_grid_constrained);
