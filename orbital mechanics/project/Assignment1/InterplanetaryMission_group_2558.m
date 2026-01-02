%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
proj_d = cd(1:backs(end)); addpath([proj_d '\student_functions']); 
addpath([proj_d '\lib']); addpath([proj_d '\lib' '\timeConversion']);
clear; close all; clc;

                        %% ============= %%
                        %% === LEG 1 === %%
                        %% ============= %%
%% 1. Initialisation
% --- Constants ---
steps = 100; mu_sun = astroConstants(4);

% --- Travel Window ---
travel_window.start_date = [2030, 1, 1, 0, 0, 0];
travel_window.end_date = [2032, 1, 1, 0, 0, 0];
travel_window.start_mjd2k = date2mjd2000(travel_window.start_date);
travel_window.end_mjd2k = date2mjd2000(travel_window.end_date);

travel_window.tof = travel_window.start_date - travel_window.end_date;

% --- Departure Planet ---
planet_mercury.name = "Mercury";
planet_mercury.id = 1;
planet_mercury.mu = astroConstants(10 + planet_mercury.id);

% --- Flyby Planet ---
planet_earth.name = "Earth";
planet_earth.id = 3;
planet_earth.mu = astroConstants(10 + planet_earth.id);

% --- Arrival Asteroid ---
asteroid.id = 316801;
asteroid.name = "N." + asteroid.id;

            %% === STATE 1/6: MERCURY === %%
            % R1 = RM1, V1 = VM1
            % RM1, VM1 known
% The satellite is at Mercury, with its state matching Mercury's state. The
% satellite's heliocentric position and velocity is the same as Mercury's
% heliocentric position and velocity. The satellite's planet-centric 
% (wrt Mercury) position is on Mercury's surface and its velocity is 0.

            %% === MANOUVRE 1/3: DEPARTURE === %%
% Lambert algorithm to solve for the manouvre which transfers the satellite
% from Mercury to Earth. This is done for every departure time to every
% arrival time, discarding the combinations where arrival time is before
% departure time of course. 

            %% === STATE 2/6: MERCURY DEPARTURE === %%
            %  --- START OF LEG 1 --- %
% The satellite is still at Mercury, but now its velocity is changed such
% that it will reach Earth at some point in the future (calculated for an 
% array of times). The satellite's heliocentric position is still the same
% as Mercury's but the its velocity is V2, the velocity given by Lambert's
% algorithm to reach Earth. The satellite's planet-centric position is
% still on its surface, and its velocity is technically v_p_plus, but it is
% irrelevant for any calculations. 

            %% === STATE 3/6: EARTH ARRIVAL === %%
            % R3 ~= RE3, V3 = V_E - v_inf_2_minus
            % RE3, VE3 known
% The satellite has now reached Earth and it has some excess velocity as it
% enters Earth's sphere of influence. The satellite's heliocentric
% position is technically the position of Earth with the size of Earth's
% sphere of influence subtracted, but this can be approximated as the
% satellite simply being at Earth's same point in the heliocentric frame.
% The heliocentric velocity is given as the second output from the Lambert
% algorithm carried out from Mercury, V3. This velocity will direct the
% satellite towards some random place in space, which is why the powered
% gravity assist is necessary. 

            %% === MANOUVRE 2/3: POWERED GRAVITY ASSIST === %%
% PLACEHOLDER

            %% === STATE 4/6: EARTH DEPARTURE === %%
            %  --- START OF LEG 2 --- %
            % R4 ~= RE4, V4 = VE4 + v_inf_2_plus
            % RE4 ~= RE3, VE4 ~= VE3
% The satellite is still technically in Earth's sphere of influence, but is
% exiting it. The heliocentric position is, again, roughly approximated as
% being equal to the Earth's heliocentric position, and the satellite's
% heliocentric velocity is the sum of the Earth's heliocentric velocity and
% the satellite's planet-centric velocity at the exit of Earth's sphere of
% influence. The satellite's planet-centric position is strictly at the
% boundary of the Earth's sphere of influence but is not relevant for any
% calculation whereas the satellite's planet-centric velocity is its
% velocity at the exit of Earth's sphere of influence (v_inf_2_plus).

            %% === STATE 5/6: ASTEROID ARRIVAL === %%
% PLACEHOLDER

            %% === MANOUVRE 3/3: ORBIT MATCHING === %%

            %% === STATE 6/6: ASTEROID ORBIT MATCHING === %%
            %  --- START OF LEG 3 --- %
% PLACEHOLDER

            %% === STATE 1/6: MERCURY === %%
% --- Define early/late departure/arrival times ---
leg1.early_dep_mjd2000 = travel_window.start_mjd2k;
leg1.early_arr_mjd2000 = travel_window.start_mjd2k;
leg1.late_dep_mjd2000 = travel_window.end_mjd2k;
leg1.late_arr_mjd2000 = travel_window.end_mjd2k;

% --- Initialise Mercury state array ---
RM_list = zeros(steps, 3);
VM_list = zeros(steps, 3);

% --- Fill Mercury state array ---
for i = 1:steps
    [RM_list(i, :), VM_list(i, :)] = ...
        get_planet_state(leg1.dep_times(i), planet_mercury.id, mu_sun);
end

            %% === STATE 2/6: MERCURY DEPARTURE === %%
% --- Initialise array for Mercury departure state ---
R2_list = zeros(steps, 3);
V2_list = zeros(steps, 3);

            %% === MANOUVRE 1/3: DEPARTURE === %%
% --- Initialise Earth state array ---
RE_list = zeros(steps, 3);
VE_list = zeros(steps, 3);

% --- Set up time array for first leg ---
leg1.dep_times = linspace(leg1.early_dep_mjd2000, leg1.late_dep_mjd2000, steps);
leg1.arr_times = linspace(leg1.early_arr_mjd2000, leg1.late_arr_mjd2000, steps);

% --- Fill Earth state array ---
for i = 1:steps
    [RE_list(i, :), VE_list(i, :)] = ...
        get_planet_state(leg1.arr_times(i), planet_earth.id, mu_sun);
end

% --- Conduct leg 1 grid search ---
disp("conducting grid search 1 (gravity-assist injection)"); tic
[V2_list, V3_list, leg1.dvtot_array, leg1.tof_array] = ...
    deep_space_injection(RM_list, VM_list, RE_list, ~, leg1.dep_times, leg1.arr_times, steps, 0);
disp("complete!"); toc

            %% === STATE 3/6: EARTH ARRIVAL === %%
% --- Initialise geocentric incoming/outgoing velocity arrays ---
v_inf_2_minus_list = V3_list - VE_list; % can be filled now
v_inf_2_plus_list = zeros(steps, 3); % must be filled in later

% also see V3_list from leg 1 grid search loop

            %% === STATE 4/6: EARTH DEPARTURE === %%
% --- Initialise heliocentric outgoing state array ---
R4_list = zeros(steps, 3);
% V4_list = zeros(steps, 3);

            %% === MANOUVRE 2/3: POWERED GRAVITY ASSIST === %%
% --- Initialise Asteroid state array ---
RE_list = zeros(steps, 3);
VE_list = zeros(steps, 3);

% --- Set up time array for first leg ---
leg1.dep_times = linspace(leg1.early_dep_mjd2000, leg1.late_dep_mjd2000, steps);
leg1.arr_times = linspace(leg1.early_arr_mjd2000, leg1.late_arr_mjd2000, steps);

% --- Fill Earth state array ---
for i = 1:steps
    [RE_list(i, :), VE_list(i, :)] = ...
        get_planet_state(leg1.arr_times(i), planet_earth.id, mu_sun);
end
% --- Conduct leg 2 grid search ---
disp("conducting grid search 2"); tic
[V4_list, V5_list, leg2.dvtot_array, leg2.tof_array] = ...
    direct_transfer(RE_list, VE_list, RA_list);
disp("complete!"); toc

            %% === STATE 5/6: ASTEROID ARRIVAL === %%
% ---

            %% === STATE 6/6: ASTEROID ORBIT MATCHING === %%
% ---

            %% === MANOUVRE 3/3: ORBIT MATCHING === %%
% ---

%% 2. Evaluate Δ𝑣tot for a grid of departure and arrival times covering ...
% the time windows provided for the Mercury-Earth leg.
% --- create array of departure and arrival times ---
leg1.dep_id = planet_mercury.id; leg1.arr_id = planet_earth.id;
leg1.early_dep_mjd2000 = travel_window.start_mjd2k;
leg1.early_arr_mjd2000 = travel_window.start_mjd2k;
leg1.late_dep_mjd2000 = travel_window.end_mjd2k;
leg1.late_arr_mjd2000 = travel_window.end_mjd2k;

leg1.dep_times = linspace(leg1.early_dep_mjd2000, leg1.late_dep_mjd2000, steps);
leg1.arr_times = linspace(leg1.early_arr_mjd2000, leg1.late_arr_mjd2000, steps);

% --- pre-calculate planet and asteroid states ---
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
    R1 = leg1.R_dep_list(i, :);
    V1 = leg1.V_dep_list(i, :);
    t1 = leg1.dep_times(i) * 24 * 3600;

    dv_row = NaN(1, steps);
    tof_row = NaN(1, steps);
    for j = 1:steps
        t2 = leg1.arr_times(j) * 24 * 3600;
        tof = t2 - t1;
        
        if tof > 0
            R3 = leg1.R_arr_list(j, :);
            V3_list = leg1.V_arr_list(j, :);

            [~, ~, ~, leg1.ERROR, V2_list, v_t2, ~, ~] = ...
                lambertMR(R1, R3, tof, mu_sun, 0, 0, 0, 0);
            if leg1.ERROR == 0
                dv_row(j) = norm(v_t1 - V1) + norm(v_t2 - V3_list);
                tof_row(j) = t2 - t1;
            end
        end
    end
    leg1.dvtot(i, :) = dv_row;
    leg1.tof(i, :) = tof_row / (24 * 3600);
end
disp("complete!"); toc

%% 3. Calculate the 


%% 4. Evaluate dvtot for a grid of departure and arrival times covering ...
% the time windows provided for the Earth-Asteroid leg

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
mjd2k1 = date2mjd2000([2031 7 15 9 41 49]);
mjd2k2 = date2mjd2000([2032 6 15 9 41 49]);
t1 = mjd2k1 * 24 * 3600;
t2 = mjd2k2 * 24 * 3600;
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);
AU = astroConstants(2);

% --- get planet and asteroid positions ---
[RD1, VD1] = get_planet_state(mjd2k1, planet_mercury.id, mu_sun);
[RGA1, VGA1] = get_planet_state(mjd2k1, planet_earth.id, mu_sun);
[RA1, VA1] = get_asteroid_state(mjd2k1, asteroid.id, mu_sun);

% c. propagate the departure planet orbit
y_D = [RD1, VD1];
tspan_D = linspace(t1, t2); % integration time span array
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_D, y_D, options);
R_D = Y(:, 1:3) ./ AU;

% d. propagate the gravity-assist planet orbit
y_GA = [RGA1, VGA1];
tspan_GA = linspace(t1, t2); % integration time span array
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_GA, y_GA, options);
R_GA = Y(:, 1:3) ./ AU;

% e. propagate the arrival asteroid orbit
y_A = [RA1, VA1];
tspan_A = linspace(t1, t2); % integration time span array
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_A, y_A, options);
R_A = Y(:, 1:3) ./ AU;

% --- plot ---
figure("Name", "Orbit Plot"); hold on;

% a. mercury-earth leg

% b. earth-asteroid leg

% c. departure planet orbit
plot3(R_D(:, 1), R_D(:, 2), R_D(:, 3), "r"); % during transfer

% d. gravity-assist planet orbit
plot3(R_GA(:, 1), R_GA(:, 2), R_GA(:, 3), "g"); % during transfer

% e. asteroid orbit
plot3(R_A(:, 1), R_A(:, 2), R_A(:, 3), "b"); % during transfer

% f. departure planet boundary positions
scatter3(R_D(1, 1), R_D(1, 2), R_D(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "r");
scatter3(R_D(end, 1), R_D(end, 2), R_D(end, 3), "filled", "MarkerFaceColor", "r");

% g. gravity-assist planet boundary positions
scatter3(R_GA(1, 1), R_GA(1, 2), R_GA(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "g");
scatter3(R_GA(end, 1), R_GA(end, 2), R_GA(end, 3), "filled", "MarkerFaceColor", "g");

% h. asteroid boundary positions
scatter3(R_A(1, 1), R_A(1, 2), R_A(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "b");
scatter3(R_A(end, 1), R_A(end, 2), R_A(end, 3), "filled", "MarkerFaceColor", "b");

% i. sun position
scatter3(0, 0, 0, "filled", "MarkerFaceColor", "y");

% --- plot properties ---
xlabel("X [AU]"); ylabel("Y [AU]"); zlabel("Z [AU]");
title("Two-body problem orbit");
legend( ...
    "", ...
    "", ...
    "", ...
    "departure planet at departure", ...
    "departure planet at arrival", ...
    "ga planet at departure", ...
    "ga planet at arrival", ...
    "asteroid at departure", ...
    "asteroid at arrival", ...
    "");
axis equal; grid on; view(3);
hold off;

%% animated plot
% 1. Define orbital data
t = linspace(t1, t2);
x = R_D(:, 1); x2 = R_GA(:, 1); x3 = R_A(:, 1);
y = R_D(:, 2); y2 = R_GA(:, 2); y3 = R_A(:, 2);
z = R_D(:, 3); z2 = R_GA(:, 3); z3 = R_A(:, 3);

% 2. Initialise the plot
pause(2)
figure("Name", "Animated Orbit Plot"); hold on; grid on; view(3); axis equal;
xlim([-max([max(abs(x)), max(abs(x2)), max(abs(x3))]), max([max(abs(x)), max(abs(x2)), max(abs(x3))])])
ylim([-max([max(abs(y)), max(abs(y2)), max(abs(y3))]), max([max(abs(y)), max(abs(y2)), max(abs(y3))])])
zlim([-max([max(abs(z)), max(abs(z2)), max(abs(z3))]), max([max(abs(z)), max(abs(z2)), max(abs(z3))])])
title("Simultaneous Multi-Orbit Animation");
scatter3(0, 0, 0, "filled", "MarkerFaceColor", "y");

% 3. Create animated lines
h1 = animatedline("Color", "r", "LineWidth", 1.5, "MaximumNumPoints", inf);
h2 = animatedline("Color", "g", "LineWidth", 1.5, "MaximumNumPoints", inf);
h3 = animatedline("Color", "b", "LineWidth", 1.5, "MaximumNumPoints", inf);

% 4. Add markers for the "heads" of the comets
head1 = plot3(x(1), y(1), z(1), "ro", "MarkerFaceColor", "r");
head2 = plot3(x2(1), y2(1), z2(1), "go", "MarkerFaceColor", "g");
head3 = plot3(x3(1), y3(1), z3(1), "bo", "MarkerFaceColor", "b");

% 5. Animation loop
for i = 1:length(t)
    % Update the tails
    addpoints(h1, x(i), y(i), z(i));
    addpoints(h2, x2(i), y2(i), z2(i));
    addpoints(h3, x3(i), y3(i), z3(i));
    
    % Update the heads
    set(head1, "XData", x(i), "YData", y(i), "ZData", z(i));
    set(head2, "XData", x2(i), "YData", y2(i), "ZData", z2(i));
    set(head3, "XData", x3(i), "YData", y3(i), "ZData", z3(i));
    
    drawnow limitrate; pause(0.05); 
end
