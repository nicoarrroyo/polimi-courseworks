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
steps = 100;

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
travel_window_close_date = [2032, 1, 1, 0, 0, 0];
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
    deep_space_injection(RM_list, VM_list, RE_list, ...
    dep_times_M, arr_times_E, steps, 0, 50);

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
    deep_space_injection(RE_list, VE_list, RA_list, ...
    dep_times_E, arr_times_A, steps, 0, 50);

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


%             %% === STATE 3/6: EARTH ARRIVAL === %%
% % --- Initialise geocentric incoming/outgoing velocity arrays ---
% v_inf_2_minus_list = V3_list - VE_list; % can be filled now
% % v_inf_2_plus_list = zeros(steps, 3); % must be filled in later
% 
% % also see V3_list from leg 1 grid search loop
% 
%                         %% ============= %%
%                         %% === LEG 2 === %%
%                         %% ============= %%
%             %% === STATE 4/6: EARTH DEPARTURE === %%
% % --- Initialise heliocentric outgoing state array ---
% R4_list = zeros(steps, steps, 3);
% V4_list = zeros(steps, steps, 3);
% 
% rp_list = zeros(steps, steps, 1);
% rp_list_invalid = [];
% rp_list_broken = [];
% 
% turn_angle = zeros(steps, steps, 1);
% 
%             %% === MANOUVRE 2/3: POWERED GRAVITY ASSIST === %%
% % --- Initialise Asteroid state array ---
% RA_list = zeros(steps, 3);
% VA_list = zeros(steps, 3);
% 
% % --- Set up time array for second leg ---
% leg2.dep_times = linspace(leg1.early_dep_mjd2000, leg1.late_dep_mjd2000, steps);
% leg2.arr_times = linspace(leg1.early_arr_mjd2000, leg1.late_arr_mjd2000, steps);
% 
% % --- Fill Asteroid state array ---
% for i = 1:steps
%     [RA_list(i, :), VA_list(i, :)] = ...
%         get_asteroid_state(leg2.arr_times(i), asteroid.id);
% end
% 
% % --- Conduct leg 2 grid search ---
% disp("conducting grid search 2 (full lambert transfer) (nested)"); tic
% for i = 1:steps
%     [V4_list(i, :, :), V5_list(i, :, :), leg2.dvtot_array, leg2.tof_array] = ...
%         deep_space_injection(RE_list, VE_list, RA_list, VA_list, leg2.dep_times, leg2.arr_times, steps, 1);
% end
% disp("complete!"); toc
% 
% % --- Calculate necessary perigee burn to match outgoing arc ---
% h_atm = 500; % height of earth atmosphere from sea-level [km]
% rp_crit = r_E + h_atm; % critical fly-by pericentre radius
% 
% disp("Calculating pericente radii for possible gravity-assist manouvres")
% for i = 1:steps
%     v_inf_2_plus_list = V4_list(i, :, :) - VE_list; % geocentric outgoing velocity
%     turn_angle(i) = acos(...
%         dot(v_inf_2_minus_list(i, :), v_inf_2_plus_list(i, :)) / ...
%         (norm(v_inf_2_minus_list(i, :)) * norm(v_inf_2_plus_list(i, :))) ...
%         );
%     if isnan(turn_angle(i))
%        continue
%     end
%     eq = @(rp) turn_angle(i) - ...
%         asin(1 / (1 + (rp * norm(v_inf_2_plus_list(i, :))^2) / mu_E)) - ...
%         asin(1 / (1 + (rp * norm(v_inf_2_minus_list(i, :))^2) / mu_E));
%     rp = fzero(eq, r_E + h_atm);
% 
%     if rp > rp_crit
%         rp_list(i) = rp;
%     elseif rp < rp_crit
%         rp_list_invalid(i) = rp;
%     else
%         rp_list_broken(i) = rp;
%     end
% end
% disp("There were " + length(rp_list_invalid) + " manouvres at invalid pericentre radii")
% disp("There were " + length(rp_list_broken) + " impossible manouvres")
% 
% ecc_2_minus_list = 1 + (rp_list .* norm(v_inf_2_minus_list).^2 / mu_E);
% a_2_minus_list = -mu_E ./ (norm(v_inf_2_plus_list).^2);
% v_p_2_minus_list = sqrt(mu_E .* (2 ./ rp_list - 1 ./ a_2_minus_list));
% 
% ecc_2_plus_list = 1 + (rp_list .* norm(v_inf_2_plus_list).^2 / mu_E);
% a_2_plus_list = -mu_E ./ (norm(v_inf_2_plus_list).^2);
% v_p_2_plus_list = sqrt(mu_E .* (2 ./ rp_list - 1 ./ a_2_plus_list));
% 
% dv_p_list = v_p_2_plus_list - v_p_2_minus_list;
% dv_ga_list = norm(dv_p_list);
% 
%             %% === STATE 5/6: ASTEROID ARRIVAL === %%
% % ---
% 
%             %% === STATE 6/6: ASTEROID ORBIT MATCHING === %%
% % ---
% 
%             %% === MANOUVRE 3/3: ORBIT MATCHING === %%
% % ---
% 
% 
% %% 5. Plot the transfer trajectory for this mission
% % --- set up times and options ---
% % mjd2k1 = date2mjd2000([2031 7 15 9 41 49]);
% % mjd2k2 = date2mjd2000([2032 6 15 9 41 49]);
% mjd2k1 = date2mjd2000([2031 7 15 9 41 49]);
% mjd2k2 = date2mjd2000([2031 10 15 9 41 49]);
% t1 = mjd2k1 * 24 * 3600;
% t2 = mjd2k2 * 24 * 3600;
% options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);
% 
% % --- get planet and asteroid positions ---
% [RD1, VD1] = get_planet_state(mjd2k1, planet_mercury.id);
% [RGA1, VGA1] = get_planet_state(mjd2k1, planet_earth.id);
% [RA1, VA1] = get_asteroid_state(mjd2k1, asteroid.id);
% 
% % c. propagate the departure planet orbit
% y_D = [RD1, VD1];
% tspan_D = linspace(t1, t2); % integration time span array
% [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_D, y_D, options);
% R_D = Y(:, 1:3) ./ AU;
% 
% % d. propagate the gravity-assist planet orbit
% y_GA = [RGA1, VGA1];
% tspan_GA = linspace(t1, t2); % integration time span array
% [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_GA, y_GA, options);
% R_GA = Y(:, 1:3) ./ AU;
% 
% % e. propagate the arrival asteroid orbit
% y_A = [RA1, VA1];
% tspan_A = linspace(t1, t2); % integration time span array
% [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_A, y_A, options);
% R_A = Y(:, 1:3) ./ AU;
% 
% % --- plot ---
% figure("Name", "Orbit Plot"); hold on;
% 
% % a. mercury-earth leg
% 
% % b. earth-asteroid leg
% 
% % c. departure planet orbit
% plot3(R_D(:, 1), R_D(:, 2), R_D(:, 3), "r"); % during transfer
% 
% % d. gravity-assist planet orbit
% plot3(R_GA(:, 1), R_GA(:, 2), R_GA(:, 3), "g"); % during transfer
% 
% % e. asteroid orbit
% plot3(R_A(:, 1), R_A(:, 2), R_A(:, 3), "b"); % during transfer
% 
% % f. departure planet boundary positions
% scatter3(R_D(1, 1), R_D(1, 2), R_D(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "r");
% scatter3(R_D(end, 1), R_D(end, 2), R_D(end, 3), "filled", "MarkerFaceColor", "r");
% 
% % g. gravity-assist planet boundary positions
% scatter3(R_GA(1, 1), R_GA(1, 2), R_GA(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "g");
% scatter3(R_GA(end, 1), R_GA(end, 2), R_GA(end, 3), "filled", "MarkerFaceColor", "g");
% 
% % h. asteroid boundary positions
% scatter3(R_A(1, 1), R_A(1, 2), R_A(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "b");
% scatter3(R_A(end, 1), R_A(end, 2), R_A(end, 3), "filled", "MarkerFaceColor", "b");
% 
% % i. sun position
% scatter3(0, 0, 0, "filled", "MarkerFaceColor", "y");
% 
% % --- plot properties ---
% xlabel("X [AU]"); ylabel("Y [AU]"); zlabel("Z [AU]");
% title("Two-body problem orbit");
% legend( ...
%     "", ...
%     "", ...
%     "", ...
%     "departure planet at departure", ...
%     "departure planet at arrival", ...
%     "ga planet at departure", ...
%     "ga planet at arrival", ...
%     "asteroid at departure", ...
%     "asteroid at arrival", ...
%     "");
% axis equal; grid on; view(3);
% hold off;

% %% animated plot
% % 1. Define orbital data
% t = linspace(t1, t2);
% x = R_D(:, 1); x2 = R_GA(:, 1); x3 = R_A(:, 1);
% y = R_D(:, 2); y2 = R_GA(:, 2); y3 = R_A(:, 2);
% z = R_D(:, 3); z2 = R_GA(:, 3); z3 = R_A(:, 3);
% 
% % 2. Initialise the plot
% pause(2)
% figure("Name", "Animated Orbit Plot"); hold on; grid on; view(3); axis equal;
% xlim([-max([max(abs(x)), max(abs(x2)), max(abs(x3))]), max([max(abs(x)), max(abs(x2)), max(abs(x3))])])
% ylim([-max([max(abs(y)), max(abs(y2)), max(abs(y3))]), max([max(abs(y)), max(abs(y2)), max(abs(y3))])])
% zlim([-max([max(abs(z)), max(abs(z2)), max(abs(z3))]), max([max(abs(z)), max(abs(z2)), max(abs(z3))])])
% title("Simultaneous Multi-Orbit Animation");
% scatter3(0, 0, 0, "filled", "MarkerFaceColor", "y");
% 
% % 3. Create animated lines
% h1 = animatedline("Color", "r", "LineWidth", 1.5, "MaximumNumPoints", inf);
% h2 = animatedline("Color", "g", "LineWidth", 1.5, "MaximumNumPoints", inf);
% h3 = animatedline("Color", "b", "LineWidth", 1.5, "MaximumNumPoints", inf);
% 
% % 4. Add markers for the "heads" of the comets
% head1 = plot3(x(1), y(1), z(1), "ro", "MarkerFaceColor", "r");
% head2 = plot3(x2(1), y2(1), z2(1), "go", "MarkerFaceColor", "g");
% head3 = plot3(x3(1), y3(1), z3(1), "bo", "MarkerFaceColor", "b");
% pause(1)
% 
% % 5. Animation loop
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
