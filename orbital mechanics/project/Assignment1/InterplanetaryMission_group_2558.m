%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
proj_d = cd(1:backs(end)); addpath([proj_d '\lib']); 
addpath([proj_d '\lib' '\timeConversion']); clear; close all; clc;

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

%% 1. Constants
steps = 200;
dv_lim = 50; % [km s^-1] should be set as low as possible (reduces computation time)

mu_sun = astroConstants(4); % Sun Gravitational Parameter [km^3 s^-2]
AU = astroConstants(2); % Astronomical Unit [km]

planet_M_id = 1;
planet_M_name = "Mercury";
planet_M_mu = astroConstants(10 + planet_M_id);

planet_E_id = 3;
planet_E_name = "Earth";
planet_E_r = astroConstants(20 + planet_E_id);
planet_E_mu = astroConstants(10 + planet_E_id);

asteroid_id = 316801;
asteroid_name = "N." + asteroid_id;

%% 2. Initialise Arrays
% --- Time ---
travel_window_start_date = [2050, 1, 1, 0, 0, 0];
travel_window_close_date = [2052, 1, 1, 0, 0, 0];
travel_window_start_mjd2k = date2mjd2000(travel_window_start_date);
travel_window_close_mjd2k = date2mjd2000(travel_window_close_date);

time_list = linspace(travel_window_start_mjd2k, travel_window_close_mjd2k, steps);

% --- Heliocentric Position (R) and Velocity (V) ---
% Mercury
RM_list = zeros(steps, 3);
VM_list = zeros(steps, 3);
for i = 1:steps
    [RM_list(i, :), VM_list(i, :)] = ...
        get_planet_state(time_list(i), planet_M_id, mu_sun);
end

% Earth
RE_list = zeros(steps, 3);
VE_list = zeros(steps, 3);
for i = 1:steps
    [RE_list(i, :), VE_list(i, :)] = ...
        get_planet_state(time_list(i), planet_E_id, mu_sun);
end

% Asteroid
RA_list = zeros(steps, 3);
VA_list = zeros(steps, 3);
for i = 1:steps
    [RA_list(i, :), VA_list(i, :)] = ...
        get_asteroid_state(time_list(i), asteroid_id, mu_sun);
end

% Satellite
V1_list = NaN(steps, steps, 3); % Departure from Mercury
V2_list = NaN(steps, steps, 3); % Arrival at Earth
V3_list = NaN(steps, steps, 3); % Departure from Earth
V4_list = NaN(steps, steps, 3); % Arrival at Asteroid

% --- Velocity Change (dv) ---
% Leg 1: Mercury-Earth
dv_grid1 = NaN(steps, steps, 3);
dv_grid1_norm = NaN(steps, steps);

% Leg 2: Earth-Asteroid
dv_grid2 = NaN(steps, steps, 3);
dv_grid2_norm = NaN(steps, steps);

% --- Time of Flight (tof) ---
tof_grid1 = NaN(steps, steps); % Leg 1: Mercury-Earth
tof_grid2 = NaN(steps, steps); % Leg 2: Earth-Asteroid

%% Manouvre 1: Gravity Assist Injection from Mercury
% --- Conduct leg 1 grid search ---
fprintf("conducting grid search 1 (gravity-assist injection)... "); tic

[V1_list, V2_list, dv_grid1, tof_grid1] = ...
    deep_space_injection(RM_list, VM_list, RE_list, VE_list, ...
    time_list, time_list, steps, 0, dv_lim);

disp("complete!");

% --- Compute dv ---
dv_grid1_norm = vecnorm(dv_grid1, 2, 3);

% --- Analyse grid search the Mercury-Earth Leg ---
find_lowest_dv_mission(dv_grid1_norm, time_list, time_list); toc

porkchop_plot( ...
    "Mercury-Earth", ...
    dv_grid1_norm, ...
    time_list, ...
    time_list);

%% Manouvre 2.1: Lambert Arc from Earth (excludes asteroid Rendez-Vous)
% --- Conduct leg 2 grid search ---
fprintf("\nconducting grid search 2 (gravity-assist)... "); tic

[V3_list, V4_list, dv_grid2, tof_grid2] = ...
    deep_space_injection(RE_list, VE_list, RA_list, VA_list, ...
    time_list, time_list, steps, 0, dv_lim);

disp("complete!");

% --- Compute dv ---
dv_grid2_norm = vecnorm(dv_grid2, 2, 3);

% --- Analyse grid search the Earth-Asteroid Leg ---
find_lowest_dv_mission(dv_grid2_norm, time_list, time_list); toc

porkchop_plot( ...
    "Earth-Asteroid", ...
    dv_grid2_norm, ...
    time_list, ...
    time_list);

%% Manouvre 2.2: Gravity Assist at Earth
[~, valid_cols_dv1] = find(~isnan(dv_grid1_norm));
valid_cols_dv1 = unique(valid_cols_dv1, "stable");

[valid_rows_dv2, ~] = find(~isnan(dv_grid2_norm));
valid_rows_dv2 = unique(valid_rows_dv2, "stable");

common_flyby_indices = intersect(valid_cols_dv1, valid_rows_dv2);

[dv1_vals, dv1_locs] = min(dv_grid1_norm(:, common_flyby_indices), [], 1, "omitnan");
[dv2_vals, dv2_locs] = min(dv_grid2_norm(common_flyby_indices, :), [], 2, "omitnan");

num_candidates = length(common_flyby_indices);
V_minus = zeros(num_candidates, 3);
V_plus = zeros(num_candidates, 3);
V_planet = zeros(num_candidates, 3); % (earth)

for k = 1:num_candidates
    t_idx = common_flyby_indices(k);
    V_minus = reshape(V2_list(dv1_locs(k), t_idx, :), 1, 3);
    V_plus = reshape(V3_list(t_idx, dv2_locs(k), :), 1, 3);
    V_planet = reshape(VE_list(t_idx, :), 1, 3);
end

v_inf_minus = V_minus - V_planet;
v_inf_plus = V_plus - V_planet;
v_inf_minus_norm = vecnorm(v_inf_minus, 2, 2);
v_inf_plus_norm = vecnorm(v_inf_plus, 2, 2);

dot_prod = sum(v_inf_minus .* v_inf_plus, 2);
delta_crit = acos(dot_prod ./ (v_inf_minus_norm .* v_inf_plus_norm));

r_crit = planet_E_r + 500; % because h_atm = 500km

term_minus  = 1 ./ (1 + (r_crit .* v_inf_minus_norm.^2) ./ planet_E_mu);
term_plus = 1 ./ (1 + (r_crit .* v_inf_plus_norm.^2) ./ planet_E_mu);

delta_max = asin(term_minus) + asin(term_plus);

is_feasible = delta_crit <= (delta_max + 10^-6); % 10^-6 for fp errors

%% 3. Stitching
total_dv_vector = dv1_vals(:)' + dv2_vals(:)'; % force row vectors for element-wise addition
total_dv_vector(~is_feasible) = Inf;
[lowest_dvtot2, best_idx] = min(total_dv_vector, [], "omitnan");

optimal_M_idx = dv1_locs(best_idx); % The row in Grid 1
optimal_E_idx = common_flyby_indices(best_idx); % The specific time index
optimal_A_idx = dv2_locs(best_idx); % The col in Grid 2

%% 4. Plot the transfer trajectory for this mission
% --- set up times and options ---
mjd2k1 = time_list(optimal_M_idx);
mjd2k2 = time_list(optimal_E_idx);
mjd2k3 = time_list(optimal_A_idx);

t1 = mjd2k1 * 24 * 3600;
t2 = mjd2k2 * 24 * 3600;
t3 = mjd2k3 * 24 * 3600;
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% --- get planet and asteroid heliocentric positions and velocities ---
[RM1, VM1] = get_planet_state(mjd2k1, planet_M_id, mu_sun);
[RE1, VE1] = get_planet_state(mjd2k1, planet_E_id, mu_sun);
[RA1, VA1] = get_asteroid_state(mjd2k1, asteroid_id, mu_sun);

% --- get satellite heliocentric positions and velocities ---
RSM = RM1;
VSM = reshape(V1_list(optimal_M_idx, optimal_E_idx, :), [1, 3]); % at mercury

[RSE, ~] = get_planet_state(mjd2k2, planet_E_id, mu_sun);
VSE = reshape(V3_list(optimal_E_idx, optimal_A_idx, :), [1, 3]); % at earth (after GA)

% a. propagate the satellite's mercury-earth leg
y_dep = [RSM, VSM];
tspan_dep = linspace(t1, t2); % integration time span array
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_dep, y_dep, options);
R_dep = Y(:, 1:3) ./ AU;

% b. propagate the satellite's earth-asteroid leg
y_GA = [RSE, VSE];
tspan_GA = linspace(t2, t3); % integration time span array
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_GA, y_GA, options);
R_GA = Y(:, 1:3) ./ AU;

% c. propagate the departure planet orbit
y_M = [RM1, VM1];
tspan_M = linspace(t1, t3); % integration time span array
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_M, y_M, options);
R_M = Y(:, 1:3) ./ AU;

% d. propagate the gravity-assist planet orbit
y_E = [RE1, VE1];
tspan_E = linspace(t1, t3); % integration time span array
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_E, y_E, options);
R_E = Y(:, 1:3) ./ AU;

% e. propagate the arrival asteroid orbit
y_A = [RA1, VA1];
tspan_A = linspace(t1, t3); % integration time span array
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_A, y_A, options);
R_A = Y(:, 1:3) ./ AU;

% --- plot ---
figure("Name", "Orbit Plot"); hold on;
% a. mercury-earth leg
plot3(R_dep(:, 1), R_dep(:, 2), R_dep(:, 3), "y--"); % during transfer

% b. earth-asteroid leg
plot3(R_GA(:, 1), R_GA(:, 2), R_GA(:, 3), "y"); % during transfer

% c. planet earth orbit
plot3(R_M(:, 1), R_M(:, 2), R_M(:, 3), "r"); % during transfer

% d. planet earth orbit
plot3(R_E(:, 1), R_E(:, 2), R_E(:, 3), "g"); % during transfer

% e. asteroid orbit
plot3(R_A(:, 1), R_A(:, 2), R_A(:, 3), "b"); % during transfer

% f. departure planet boundary positions
scatter3(R_M(1, 1), R_M(1, 2), R_M(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "r");
scatter3(R_M(end, 1), R_M(end, 2), R_M(end, 3), "filled", "MarkerFaceColor", "r");

% g. gravity-assist planet boundary positions
scatter3(R_E(1, 1), R_E(1, 2), R_E(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "g");
scatter3(R_E(end, 1), R_E(end, 2), R_E(end, 3), "filled", "MarkerFaceColor", "g");

% h. asteroid boundary positions
scatter3(R_A(1, 1), R_A(1, 2), R_A(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "b");
scatter3(R_A(end, 1), R_A(end, 2), R_A(end, 3), "filled", "MarkerFaceColor", "b");

% i. sun position
scatter3(0, 0, 0, "filled", "MarkerFaceColor", "y");

% --- plot properties ---
xlabel("X [AU]"); ylabel("Y [AU]"); zlabel("Z [AU]");
title("Two-body problem orbit");
legend( ...
    "", ... a. mercury-earth leg
    "", ... b. earth-asteroid leg
    "", ... c. planet mercury orbit
    "", ... d. planet earth orbit
    "", ... e. asteroid orbit
    "departure planet at departure", ...
    "departure planet at arrival", ...
    "ga planet at departure", ...
    "ga planet at arrival", ...
    "asteroid at departure", ...
    "asteroid at arrival", ...
    "" ... i. sun position
    );
axis equal; grid on; view(3);
hold off;

%% 5. animated plot
% % Define orbital data
% t = linspace(t1, t2);
% x = R_M(:, 1); x2 = R_E(:, 1); x3 = R_A(:, 1); 
% y = R_M(:, 2); y2 = R_E(:, 2); y3 = R_A(:, 2);
% z = R_M(:, 3); z2 = R_E(:, 3); z3 = R_A(:, 3);
% 
% % Initialise the plot
% pause(2)
% figure("Name", "Animated Orbit Plot"); hold on; grid on; view(3); axis equal;
% xlim([-max([max(abs(x)), max(abs(x2)), max(abs(x3))]), max([max(abs(x)), max(abs(x2)), max(abs(x3))])])
% ylim([-max([max(abs(y)), max(abs(y2)), max(abs(y3))]), max([max(abs(y)), max(abs(y2)), max(abs(y3))])])
% zlim([-max([max(abs(z)), max(abs(z2)), max(abs(z3))]), max([max(abs(z)), max(abs(z2)), max(abs(z3))])])
% title("Simultaneous Multi-Orbit Animation");
% scatter3(0, 0, 0, "filled", "MarkerFaceColor", "y");
% 
% % Create animated lines
% h1 = animatedline("Color", "r", "LineWidth", 1.5, "MaximumNumPoints", inf);
% h2 = animatedline("Color", "g", "LineWidth", 1.5, "MaximumNumPoints", inf);
% h3 = animatedline("Color", "b", "LineWidth", 1.5, "MaximumNumPoints", inf);
% 
% % Add markers for the "heads" of the comets
% head1 = plot3(x(1), y(1), z(1), "ro", "MarkerFaceColor", "r");
% head2 = plot3(x2(1), y2(1), z2(1), "go", "MarkerFaceColor", "g");
% head3 = plot3(x3(1), y3(1), z3(1), "bo", "MarkerFaceColor", "b");
% pause(1)
% 
% % Animation loop
% for i = 1:length(t)
%     % Update the tails
%     addpoints(h1, x(i), y(i), z(i));
%     addpoints(h2, x2(i), y2(i), z2(i));
%     addpoints(h3, x3(i), y3(i), z3(i));
% 
%     % Update the heads
%     set(head1, "XData", x(i), "YData", y(i), "ZData", z(i));
%     set(head2, "XData", x2(i), "YData", y2(i), "ZData", z2(i));
%     set(head3, "XData", x3(i), "YData", y3(i), "ZData", z3(i));
% 
%     drawnow limitrate; pause(0.05); 
% end
