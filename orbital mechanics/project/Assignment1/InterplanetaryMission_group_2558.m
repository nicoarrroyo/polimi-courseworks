%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
proj_d = cd(1:backs(end)); addpath([proj_d '\student_functions']); 
addpath([proj_d '\lib']); addpath([proj_d '\lib' '\timeConversion']);
clear; close all; clc;

%% ============= %%
%% === LEG 1 === %%
%% ============= %%
%% --- Mission Requirements ---
% --- Constants ---
steps = 100; mu_sun = astroConstants(4);

% --- Travel Window ---
travel_window.start_date = [2030, 1, 1, 0, 0, 0];
travel_window.end_date = [2032, 1, 1, 0, 0, 0];
travel_window.start_mjd2k = date2mjd2000(travel_window.start_date);
travel_window.end_mjd2k = date2mjd2000(travel_window.end_date);

travel_window.tof = travel_window.start_date - travel_window.end_date;

% --- Departure Planet ---
planet_dep.name = "Mercury";
planet_dep.id = 1;
planet_dep.mu = astroConstants(10 + planet_dep.id);

% --- Flyby Planet ---
planet_fb.name = "Earth";
planet_fb.id = 3;
planet_fb.mu = astroConstants(10 + planet_fb.id);

% --- Arrival Asteroid ---
asteroid_arr.id = 316801;
asteroid_arr.name = "N." + asteroid_arr.id;
[...
    asteroid_arr.kep1, ...
    asteroid_arr.mass, ...
    asteroid_arr.M1...
    ] = ephAsteroids(travel_window.start_mjd2k, asteroid_arr.id);
[...
    asteroid_arr.kep2, ...
    ~, ...
    asteroid_arr.M2...
    ] = ephAsteroids(travel_window.end_mjd2k, asteroid_arr.id);

%% 2. Evaluate Δ𝑣tot for a grid of departure and arrival times covering ...
% the time windows provided for the Mercury-Earth leg.
% --- create array of departure and arrival times ---
leg1.dep_id = planet_dep.id; leg1.arr_id = planet_fb.id;
leg1.early_dept_mjd2000 = travel_window.start_mjd2k;
leg1.early_arr_mjd2000 = travel_window.start_mjd2k;
leg1.late_dept_mjd2000 = travel_window.end_mjd2k;
leg1.late_arr_mjd2000 = travel_window.end_mjd2k;

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
leg1.v_inf_minus = NaN(steps, steps);
leg1.v_inf_plus = NaN(steps, steps);
disp("conducting grid search 1"); tic
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
leg1.dv_min = 15; leg1.dv_max = 100;
leg1.v_levels = leg1.dv_min : 1 : leg1.dv_max;
contour(leg1.dep_times, leg1.arr_times, leg1.dvtot', leg1.v_levels, "ShowText", "off");
clim([leg1.dv_min, leg1.dv_max]);

% --- plot the constant time of flight lines ---
% leg1.tof_levels = 50 : 50 : 200;
% contour(leg1.dep_times, leg1.arr_times, leg1.tof, leg1.tof_levels, "ShowText", "on");

% --- plot ---
colorbar;
xlabel("Departure Time (MJD2000)");
ylabel("Arrival Time (MJD2000)");
title("Porkchop Plot: Δv_{tot} for Mercury-Earth Leg");
hold off;

%% 4. Find the cheapest mission in terms of Δv
[leg1.min_val, min_idx_linear] = min(leg1.dvtot(:));
[row_idx, col_idx] = ind2sub(size(leg1.dvtot), min_idx_linear);

leg1.mjd_opt_dep_grid = leg1.dep_times(row_idx);
leg1.date_d = mjd20002date(leg1.mjd_opt_dep_grid);

leg1.mjd_opt_arr_grid = leg1.arr_times(col_idx);
leg1.date_a = mjd20002date(leg1.mjd_opt_arr_grid);

leg1.tof_opt_grid = leg1.mjd_opt_arr_grid - leg1.mjd_opt_dep_grid;

% --- results output ---
fprintf("\n=== GRID SEARCH RESULTS ===\n");
fprintf("Min Delta V: %.4f km s^-1\n", leg1.min_val);
fprintf("Departure:   MJD2000 %.2f\n", leg1.mjd_opt_dep_grid);
fprintf("             Date    %.0f %.0f %.0f %.0f %.0f %.0f\n", leg1.date_d);
fprintf("Arrival:     MJD2000 %.2f\n", leg1.mjd_opt_arr_grid);
fprintf("             Date    %.0f %.0f %.0f %.0f %.0f %.0f\n", leg1.date_a);
fprintf("Total ToF (days):    %.2f\n", leg1.tof_opt_grid);

%% ============= %%
%% === LEG 2 === %%
%% ============= %%
%% 2. Evaluate Δ𝑣tot for a grid of departure and arrival times covering ...
% the time windows provided for the Earth-asteroid leg.
% --- create array of departure and arrival times ---
leg2.dep_id = leg1.arr_id; leg2.arr_id = 316801;
leg2.early_dept_mjd2000 = leg1.early_arr_mjd2000;
leg2.early_arr_mjd2000 = leg1.early_arr_mjd2000; % unrealistically short tof
leg2.late_dept_mjd2000 = leg1.late_arr_mjd2000; % unrealistically short tof
leg2.late_arr_mjd2000 = leg1.late_arr_mjd2000;

leg2.dep_times = linspace(leg2.early_dept_mjd2000, leg2.late_dept_mjd2000, steps);
leg2.arr_times = linspace(leg2.early_arr_mjd2000, leg2.late_arr_mjd2000, steps);

% --- pre-calculate planet states ---
leg2.R_dep_list = leg1.R_arr_list;
leg2.V_dep_list = leg1.V_arr_list;
leg2.R_arr_list = zeros(steps, 3);
leg2.V_arr_list = zeros(steps, 3);

for i = 1:steps
    [leg2.R_arr_list(i, :), leg2.V_arr_list(i, :)] = ...
        get_asteroid_state(leg2.dep_times(i), leg2.arr_id, mu_sun);
end
[asteroid_arr.kep1, ~, asteroid_arr.M1] = ephAsteroids(leg2.dep_times(i), leg2.arr_id);

% --- conduct grid search ---
leg2.dvtot = NaN(steps, steps);
leg2.tof = NaN(steps, steps);
disp("conducting grid search 2"); tic
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
leg2.dv_min = 15; leg2.dv_max = 60;
leg2.v_levels = leg2.dv_min : 1 : leg2.dv_max;
contour(leg2.dep_times, leg2.arr_times, leg2.dvtot', leg2.v_levels, "ShowText", "off");
clim([leg2.dv_min, leg2.dv_max]);

% --- plot the constant time of flight lines ---
% leg2.tof_levels = 200 : 100 : 500;
% contour(leg2.dep_times, leg2.arr_times, leg2.tof, leg2.tof_levels, "ShowText", "on");

% --- plot ---
colorbar;
xlabel("Departure Time (MJD2000)");
ylabel("Arrival Time (MJD2000)");
title("Porkchop Plot: Δv_{tot} for Mercury-Earth Leg");
hold off;

%% 4. Find the cheapest mission in terms of Δv
[leg2.min_val, min_idx_linear] = min(leg2.dvtot(:));
[row_idx, col_idx] = ind2sub(size(leg2.dvtot), min_idx_linear);

leg2.mjd_opt_dep_grid = leg2.dep_times(row_idx);
leg2.date_d = mjd20002date(leg2.mjd_opt_dep_grid);

leg2.mjd_opt_arr_grid = leg2.arr_times(col_idx);
leg2.date_a = mjd20002date(leg2.mjd_opt_arr_grid);

leg2.tof_opt_grid = leg2.mjd_opt_arr_grid - leg2.mjd_opt_dep_grid;

% --- results output ---
fprintf("\n=== GRID SEARCH RESULTS ===\n");
fprintf("Min Delta V: %.4f km s^-1\n", leg2.min_val);
fprintf("Departure:   MJD2000 %.2f\n", leg2.mjd_opt_dep_grid);
fprintf("             Date    %.0f %.0f %.0f %.0f %.0f %.0f\n", leg2.date_d);
fprintf("Arrival:     MJD2000 %.2f\n", leg2.mjd_opt_arr_grid);
fprintf("             Date    %.0f %.0f %.0f %.0f %.0f %.0f\n", leg2.date_a);
fprintf("Total ToF (days):    %.2f\n", leg2.tof_opt_grid);


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
