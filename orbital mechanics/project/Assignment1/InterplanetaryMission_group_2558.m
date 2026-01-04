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
steps = 1000;
dv_lim = 40; % [km s^-1] should be set as low as possible (reduces computation time)

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
travel_window_start_date = [2030, 1, 1, 0, 0, 0];
travel_window_close_date = [2060, 1, 1, 0, 0, 0];
travel_window_start_mjd2k = date2mjd2000(travel_window_start_date);
travel_window_close_mjd2k = date2mjd2000(travel_window_close_date);

dep_times_M = linspace(travel_window_start_mjd2k, travel_window_close_mjd2k, steps);
arr_times_E = linspace(travel_window_start_mjd2k, travel_window_close_mjd2k, steps);
dep_times_E = linspace(travel_window_start_mjd2k, travel_window_close_mjd2k, steps);
arr_times_A = linspace(travel_window_start_mjd2k, travel_window_close_mjd2k, steps);

% --- Heliocentric Position (R) and Velocity (V) ---
% Mercury
RM_list = zeros(steps, 3);
VM_list = zeros(steps, 3);
for i = 1:steps
    [RM_list(i, :), VM_list(i, :)] = ...
        get_planet_state(dep_times_M(i), planet_M_id, mu_sun);
end

% Earth
RE_list = zeros(steps, 3);
VE_list = zeros(steps, 3);
for i = 1:steps
    [RE_list(i, :), VE_list(i, :)] = ...
        get_planet_state(dep_times_E(i), planet_E_id, mu_sun);
end

% Asteroid
RA_list = zeros(steps, 3);
VA_list = zeros(steps, 3);
for i = 1:steps
    [RA_list(i, :), VA_list(i, :)] = ...
        get_asteroid_state(arr_times_A(i), asteroid_id, mu_sun);
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
    dep_times_M, arr_times_E, steps, 0, dv_lim);

disp("complete!");

% --- Compute dv ---
dv_grid1_norm = vecnorm(dv_grid1, 2, 3);

% --- Analyse grid search the Mercury-Earth Leg ---
find_lowest_dv_mission(dv_grid1_norm, dep_times_M, arr_times_E); toc

porkchop_plot( ...
    "Mercury-Earth", ...
    dv_grid1_norm, ...
    dep_times_M, ...
    arr_times_E);

%% Manouvre 2.1: Lambert Arc from Earth
% --- Incoming Geocentric Velocity ---
v_minus_list = V2_list - reshape(VE_list, [size(V2_list, 1), 1, size(V2_list, 3)]);

% --- Conduct leg 2 grid search ---
fprintf("\nconducting grid search 2 (gravity-assist)... "); tic

[V3_list, V4_list, dv_grid2, tof_grid2] = ...
    deep_space_injection(RE_list, VE_list, RA_list, VA_list, ...
    dep_times_E, arr_times_A, steps, 1, dv_lim);

disp("complete!");

% --- Compute dv ---
dv_grid2_norm = vecnorm(dv_grid2, 2, 3);

% --- Outgoing Geocentric Velocity ---
v_plus_list = V3_list - reshape(VE_list, [size(V3_list, 1), 1, size(V3_list, 3)]);

% --- Analyse grid search the Earth-Asteroid Leg ---
find_lowest_dv_mission(dv_grid2_norm, dep_times_E, arr_times_A); toc

porkchop_plot( ...
    "Earth-Asteroid", ...
    dv_grid2_norm, ...
    dep_times_E, ...
    arr_times_A);

%% === STITCHING & OPTIMISATION === %%
% This section replaces the manual "Manouvre 2.2" and "Manouvre 3" checks.
% It uses the grid search results as a starting point to find the exact 
% continuous trajectory.

fprintf("\n=== INITIALISING GLOBAL OPTIMISER ===\n");

% 1. Extract Best Dates from Grid Search Results
% (We re-find the minimum indices here to capture the specific dates)

% --- Leg 1 Best Guess ---
[min_val1, idx1] = min(dv_grid1_norm(:));
[r1, c1] = ind2sub(size(dv_grid1_norm), idx1);
t_dep_M_guess = dep_times_M(c1);
t_arr_E_guess = arr_times_E(r1);

% --- Leg 2 Best Guess ---
[min_val2, idx2] = min(dv_grid2_norm(:));
[r2, c2] = ind2sub(size(dv_grid2_norm), idx2);
t_dep_E_guess = dep_times_E(c2);
t_arr_A_guess = arr_times_A(r2);

% 2. Create the "Stitched" Initial Guess (x0)
% We average the arrival at Earth (Leg 1) and Departure from Earth (Leg 2)
% to find a single "Flyby Date".
t_flyby_guess = (t_arr_E_guess + t_dep_E_guess) / 2;

x0 = [t_dep_M_guess, t_flyby_guess, t_arr_A_guess];

fprintf("Initial Guess Dates [MJD2000]:\n");
fprintf("Dep Mercury: %.2f\n   Flyby Earth: %.2f\n   Arr Asteroid: %.2f\n", x0);

% 3. Run the Optimiser
% We define a handle to the cost function below, passing all fixed constants.
cost_func = @(x) mission_cost(x, planet_M_id, planet_E_id, asteroid_id, mu_sun, planet_E_r, planet_E_mu);

% Options for fminsearch (display iterations so you can see progress)
options = optimset('Display', 'off', 'TolX', 1e-4, 'MaxFunEvals', 1500);

fprintf("\nRunning fminsearch (Nelder-Mead)...\n");
tic;
[x_opt, min_dv_total] = fminsearch(cost_func, x0, options);
toc;

%% === FINAL OUTPUT === %%

% Convert optimised MJD back to Dates for readability
date_dep = mjd20002date(x_opt(1));
date_fb  = mjd20002date(x_opt(2));
date_arr = mjd20002date(x_opt(3));

fprintf("\n============================================\n");
fprintf("      MISSION OPTIMISATION COMPLETE       \n");
fprintf("============================================\n");
fprintf("Total Mission Delta-V: %.4f km/s\n", min_dv_total);
fprintf("--------------------------------------------\n");
fprintf("Departure (Mercury):   %s (MJD %.2f)\n", datestr(date_dep), x_opt(1));
fprintf("Flyby (Earth):         %s (MJD %.2f)\n", datestr(date_fb), x_opt(2));
fprintf("Arrival (Asteroid):    %s (MJD %.2f)\n", datestr(date_arr), x_opt(3));
fprintf("============================================\n");

% -----------------------------------------------------------
% LOCAL FUNCTION: Mission Cost Calculator
% -----------------------------------------------------------
function J = mission_cost(x, id_dep, id_fb, id_arr, mu, R_planet_fb, mu_planet_fb)
    % Unpack decision vector
    t1 = x(1); % Departure
    t2 = x(2); % Flyby
    t3 = x(3); % Arrival
    
    % Ensure time moves forward (penalty if t2 < t1 or t3 < t2)
    if t2 <= t1 || t3 <= t2
        J = 1e5; return; 
    end
    
    % --- 1. Get Ephemerides ---
    [r1, v1_p] = get_planet_state(t1, id_dep, mu);
    [r2, v2_p] = get_planet_state(t2, id_fb, mu);
    [r3, v3_a] = get_asteroid_state(t3, id_arr, mu);
    
    % --- 2. Leg 1: Mercury -> Earth ---
    tof1 = (t2 - t1) * 86400;
    % Using lambertMR (assuming it's in your path)
    [~, ~, ~, err1, v1_dep, v1_arr, ~, ~] = lambertMR(r1, r2, tof1, mu, 0, 0, 0, 0);
    
    % --- 3. Leg 2: Earth -> Asteroid ---
    tof2 = (t3 - t2) * 86400;
    [~, ~, ~, err2, v2_dep, v2_arr, ~, ~] = lambertMR(r2, r3, tof2, mu, 0, 0, 0, 0);
    
    if err1 ~= 0 || err2 ~= 0
        J = 1e5; return; % Penalty for solver failure
    end
    
    % --- 4. Compute Costs ---
    
    % A. Launch Cost (Relative to Mercury)
    dv_launch = norm(v1_dep - v1_p);
    
    % B. Arrival Cost (Relative to Asteroid)
    dv_arr = norm(v2_arr - v3_a);
    
    % C. Powered Flyby Cost (At Earth)
    v_inf_in  = v1_arr - v2_p; % Velocity relative to Earth (Incoming)
    v_inf_out = v2_dep - v2_p; % Velocity relative to Earth (Outgoing)
    
    % Check Turn Angle Strategy
    % For the optimiser, we use a "Powered Flyby Vector Difference" approximation.
    % This drives the optimiser to align the vectors as much as possible.
    % (Refining this to a full pericentre burn model can be done, but 
    % vector difference is smoother for fminsearch convergence).
    
    dv_flyby = norm(v_inf_out - v_inf_in);
    
    % Total Cost Function
    J = dv_launch + dv_flyby + dv_arr;
end
