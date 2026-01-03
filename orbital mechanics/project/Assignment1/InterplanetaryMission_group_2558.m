%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
proj_d = cd(1:backs(end)); addpath([proj_d '\student_functions']); 
addpath([proj_d '\lib']); addpath([proj_d '\lib' '\timeConversion']);
clear; close all; clc;

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
dv_lim = 20; % [km s^-1] should be set as low as possible (reduces computation time)

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
travel_window_start_date = [2051, 1, 1, 0, 0, 0];
travel_window_close_date = [2052, 1, 1, 0, 0, 0];
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
        get_planet_state(dep_times_M(i), planet_M_id);
end

% Earth
RE_list = zeros(steps, 3);
VE_list = zeros(steps, 3);
for i = 1:steps
    [RE_list(i, :), VE_list(i, :)] = ...
        get_planet_state(dep_times_E(i), planet_E_id);
end

% Asteroid
RA_list = zeros(steps, 3);
VA_list = zeros(steps, 3);
for i = 1:steps
    [RA_list(i, :), VA_list(i, :)] = ...
        get_asteroid_state(arr_times_A(i), asteroid_id);
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
disp("conducting grid search 1 (gravity-assist injection)"); tic

[V1_list, V2_list, dv_grid1, tof_grid1] = ...
    deep_space_injection(RM_list, VM_list, RE_list, VE_list, ...
    dep_times_M, arr_times_E, steps, 0, dv_lim);

disp("complete!"); toc

% --- Compute dv ---
dv_grid1_norm = vecnorm(dv_grid1, 2, 3);

% --- Analyse grid search the Mercury-Earth Leg ---
porkchop_plot( ...
    "Mercury-Earth", ...
    dv_grid1_norm, ...
    dep_times_M, ...
    arr_times_E);

find_lowest_dv_mission(dv_grid1_norm, dep_times_M, arr_times_E);

%% Manouvre 2.1: Lambert Arc from Earth
% --- Incoming Geocentric Velocity ---
v_minus_list = V2_list - reshape(VE_list, [size(V2_list, 1), 1, size(V2_list, 3)]);

% --- Conduct leg 2 grid search ---
disp("conducting grid search 2 (gravity-assist)"); tic

[V3_list, V4_list, dv_grid2, tof_grid2] = ...
    deep_space_injection(RE_list, VE_list, RA_list, VA_list, ...
    dep_times_E, arr_times_A, steps, 1, dv_lim);

disp("complete"); toc

% --- Compute dv ---
dv_grid2_norm = vecnorm(dv_grid2, 2, 3);

% --- Outgoing Geocentric Velocity ---
v_plus_list = V3_list - reshape(VE_list, [size(V3_list, 1), 1, size(V3_list, 3)]);

% --- Analyse grid search the Earth-Asteroid Leg ---
porkchop_plot( ...
    "Earth-Asteroid", ...
    dv_grid2_norm, ...
    dep_times_E, ...
    arr_times_A);

find_lowest_dv_mission(dv_grid2_norm, dep_times_E, arr_times_A);

%% Manouvre 2.2: Gravity Assist at Earth
v_minus_norm = vecnorm(v_minus_list, 2, 3);
v_plus_norm = vecnorm(v_plus_list, 2, 3);

dot_prod = sum(v_minus_list .* v_plus_list, 3);
turn_angle = acos(dot_prod ./ (v_minus_norm .* v_plus_norm));

valid_indices = find(~isnan(turn_angle));

h_atm = 500; % Earth atmosphere altitude [km]
rp_crit = planet_E_r + h_atm; % lowest allowed fly-by radius
options = optimset("Display", "off"); % suppress warnings and iterations

rp_list = NaN(size(turn_angle));
rp_list_low = NaN(size(turn_angle));
rp_list_broken = NaN(size(turn_angle));

disp("conducting non-linear pericentre radius search"); tic;
for k = 1:length(valid_indices)
    idx = valid_indices(k);
    
    current_turn = turn_angle(idx);
    v_p_n = v_plus_norm(idx);
    v_m_n = v_minus_norm(idx);
    
    eq = @(rp) current_turn - ...
        asin(1 / (1 + (rp * v_p_n^2) / planet_E_mu)) - ...
        asin(1 / (1 + (rp * v_m_n^2) / planet_E_mu));
    
    try
        temp_rp = fzero(eq, rp_crit, options);
    catch % non-convergence case
        continue;
    end
    
    if temp_rp >= rp_crit
        rp_list(idx) = temp_rp;
    elseif temp_rp < rp_crit
        rp_list_low(idx) = temp_rp;
    else
        rp_list_broken(idx) = temp_rp;
    end
end
disp("complete!"); toc;

%% Manouvre 3: Rendez-Vous at Asteroid

