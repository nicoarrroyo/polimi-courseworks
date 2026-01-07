%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
proj_d = cd(1:backs(end)); addpath([proj_d '\lib']); 
addpath([proj_d '\lib' '\timeConversion']); clear; close all; clc;

%% 1. Constants
steps = 200;
dv_lim = 30;

mu_sun = astroConstants(4); % Sun Gravitational Parameter [km^3 s^-2]
AU = astroConstants(2); % Astronomical Unit [km]

planet_M_id = 1;
planet_M_mu = astroConstants(10 + planet_M_id);

planet_E_id = 3;
planet_E_r = astroConstants(20 + planet_E_id);
planet_E_mu = astroConstants(10 + planet_E_id);

asteroid_id = 316801;
asteroid_name = "N." + asteroid_id;

%% 2. Initialise Arrays
% --- Time ---
travel_window_start_date = [2030, 1, 1, 0, 0, 0];
travel_window_close_date = [2035, 1, 1, 0, 0, 0];
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
V1_grid = NaN(steps, steps, 3); % Departure from Mercury
V2_grid = NaN(steps, steps, 3); % Arrival at Earth
V3_grid = NaN(steps, steps, 3); % Departure from Earth
V4_grid = NaN(steps, steps, 3); % Arrival at Asteroid

% --- Velocity Change (dv) ---
% Manouvre 1: Mercury-Earth
dv_grid1 = Inf(steps, steps, 3);
dv_grid1_norm = Inf(steps, steps);

% Manouvre 2: Earth-Asteroid Gravity Assist
% variables defined when needed. no need for arrays

% Manouvre 3: Asteroid Rendez-Vous
dv_grid3 = Inf(steps, steps, 3);
dv_grid3_norm = Inf(steps, steps);

% --- Time of Flight (tof) ---
tof_grid1 = NaN(steps, steps); % Leg 1: Mercury-Earth
tof_grid3 = NaN(steps, steps); % Leg 2: Earth-Asteroid

%% Manouvre 1: Gravity Assist Injection from Mercury
% --- Conduct leg 1 grid search ---
fprintf("conducting grid search 1 (gravity-assist injection)... "); tic        

for i = 1:steps
    t1 = time_list(i) * 24 * 3600;
    R1 = RM_list(i, :);
    
    for j = 1:steps
        t2 = time_list(j) * 24 * 3600;
        tof = t2 - t1;

        % check for arrival being after departure
        if tof <= (10*24*3600) % minimum time for tof
            continue
        end

        R2 = RE_list(j, :);
        [~, ~, ~, ERROR, V_dep, V_arr, ~, ~] = ...
            lambertMR(R1, R2, tof, mu_sun, 0, 0, 0, 0);
        
        % check for valid lambert arc
        if ERROR ~= 0
            continue
        end

        dv = V_dep - VM_list(i, :); % departure from mercury

        % check for reasonable delta-v
        if norm(dv) > dv_lim
            continue
        end

        V1_grid(i, j, :) = V_dep;
        V2_grid(i, j, :) = V_arr;
        dv_grid1(i, j, :) = dv;
        tof_grid1(i, j) = tof;
    end
end

% --- Compute dv ---
dv_grid1_norm = vecnorm(dv_grid1, 2, 3);

% --- Find valid options ---
% rows = mercury departure, columns = earth arrival
[dv_grid1_valid_rows, dv_grid1_valid_cols] = find(~isnan(dv_grid1_norm));
dv_grid1_valid_rows = unique(dv_grid1_valid_rows, "stable");
dv_grid1_valid_cols = unique(dv_grid1_valid_cols, "stable");
disp("complete!");
disp("found " + length(dv_grid1_valid_rows)*length(dv_grid1_valid_cols) + " valid options out of " + (steps*steps));
toc

% --- Analyse grid search the Mercury-Earth Leg ---
porkchop_plot( ...
    "Mercury-Earth", ...
    dv_grid1_norm, ...
    time_list, ...
    time_list);

%% Manouvre 3: Lambert Arc from Earth (includes asteroid rendez-vous)
% manouvre 3 done before manouvre 2 because of problem geometry.
% --- Conduct leg 2 grid search ---
fprintf("\nconducting grid search 2 (asteroid arrival)... "); tic

for i = 1:length(dv_grid1_valid_cols)
    ii = dv_grid1_valid_cols(i);
    t1 = time_list(ii) * 24 * 3600;
    R1 = RE_list(ii, :);
    
    for j = dv_grid1_valid_cols(1):steps
        t2 = time_list(j) * 24 * 3600;
        tof = t2 - t1;

        % check for arrival being after departure
        if tof <= (10*24*3600) % minimum time for tof
            continue
        end

        R2 = RA_list(j, :);
        [~, ~, ~, ERROR, V1, V2, ~, ~] = ...
            lambertMR(R1, R2, tof, mu_sun, 0, 0, 0, 0);
        
        % check for valid lambert arc
        if ERROR ~= 0
            continue
        end

        dv = VA_list(j, :) - V2; % orbit-match at asteroid

        % check for reasonable delta-v
        if norm(dv) > dv_lim
            continue
        end

        V3_grid(ii, j, :) = V1;
        V4_grid(ii, j, :) = V2;
        dv_grid3(ii, j, :) = dv;
        tof_grid3(ii, j) = tof;
    end
end

% --- Compute dv ---
dv_grid3_norm = vecnorm(dv_grid3, 2, 3);

% --- Find valid options ---
% rows = earth departure, columns = asteroid arrival
[dv_grid3_valid_rows, dv_grid3_valid_cols] = find(~isnan(dv_grid3_norm));
dv_grid3_valid_rows = unique(dv_grid3_valid_rows, "stable");
dv_grid3_valid_cols = unique(dv_grid3_valid_cols, "stable");
disp("complete!");
disp("found " + length(dv_grid3_valid_rows)*length(dv_grid3_valid_cols) + " valid options out of " + (length(dv_grid1_valid_rows)*steps));
toc

% --- Analyse grid search the Earth-Asteroid Leg ---
porkchop_plot( ...
    "Earth-Asteroid", ...
    dv_grid3_norm, ...
    time_list, ...
    time_list);

%% Manouvre 2: Gravity Assist at Earth
fprintf("\nconducting grid search 3 (flyby)... "); tic
possible_flyby_idxs = intersect(dv_grid1_valid_cols, dv_grid3_valid_rows);
rp_crit = planet_E_r + 500;

opt_dv_tot = Inf;
for j = 1:length(possible_flyby_idxs) % for each valid "being at earth"
    jj = possible_flyby_idxs(j);
    V_planet = VE_list(jj, :);
    valid_depart = find(~isnan(V2_grid(:, jj)));
    for i = 1:length(valid_depart) % for each valid mercury departure
        ii = valid_depart(i);
        V_minus = reshape(V2_grid(ii, jj, :), 1, 3);

        v_inf_minus = V_minus - V_planet;
        v_inf_minus_norm = norm(v_inf_minus);
        e_minus_crit = 1 + rp_crit*v_inf_minus_norm^2/planet_E_mu;
        delta_minus_crit = 2*asin(1/e_minus_crit);

        this_opt_dv_tot = Inf;
        dv_launch_norm = dv_grid1_norm(ii, jj);

        for k = 1:length(dv_grid3_valid_cols) % for each valid asteroid arrival
            kk = dv_grid3_valid_cols(k);
            V_plus = reshape(V3_grid(jj, kk, :), 1, 3);

            dv_fb = V_plus - V_minus;
            dv_fb_norm = norm(dv_fb);
            if dv_fb_norm > dv_lim
                continue
            end

            dv_tot_norm = dv_launch_norm + dv_fb_norm + dv_grid3_norm(jj, kk);
            if dv_tot_norm > this_opt_dv_tot
                continue
            end

            v_inf_plus = V_plus - V_planet;
            v_inf_plus_norm = norm(v_inf_plus);

            dot_prod = dot(v_inf_minus, v_inf_plus);
            delta = acos(dot_prod / (v_inf_minus_norm * v_inf_plus_norm));

            e_plus_crit = 1 + rp_crit*v_inf_minus_norm^2/planet_E_mu;
            delta_plus_crit = 2*asin(1/e_plus_crit);
            delta_crit = (delta_minus_crit+delta_plus_crit)/2;

            if delta > delta_crit || isnan(delta)
                continue
            end

            if dv_tot_norm < norm(opt_dv_tot)
                opt_dv_launch = reshape(dv_grid1(ii, jj, :), 1, 3);
                opt_dv_fb = dv_fb;
                opt_dv_rv = reshape(dv_grid3(jj, kk, :), 1, 3);

                opt_dv_launch_norm = norm(opt_dv_launch);
                opt_dv_fb_norm = dv_fb_norm;
                opt_dv_rv_norm = norm(opt_dv_rv);

                opt_dv_tot = opt_dv_launch + opt_dv_fb + opt_dv_rv;
                opt_dv_tot_norm = opt_dv_launch_norm + opt_dv_fb_norm + opt_dv_rv_norm;

                opt_delta = delta;
                opt_v_inf_plus_norm = v_inf_plus_norm;
                opt_v_inf_minus_norm = v_inf_minus_norm;
            end
        end
    end
end

eq = @(rp) opt_delta - ...
    asin(1 / (1 + (rp * opt_v_inf_plus_norm^2) / planet_E_mu)) - ...
    asin(1 / (1 + (rp * opt_v_inf_minus_norm^2) / planet_E_mu));
rp_ans = fzero(eq, rp_crit, optimset("Display", "off"));

disp("complete!");
toc

%% Final Results Output
[optimal_M_idx, optimal_E_idx] = find(dv_grid1_norm == opt_dv_launch_norm);
[~, optimal_A_idx] = find(dv_grid3_norm == opt_dv_rv_norm);

fprintf("\nTOTAL ΔV REQUIRED: %.4f km s^-1\n", opt_dv_tot_norm);
fprintf("LEG 1 ΔV: %.4f km s^-1\n", opt_dv_launch_norm);
fprintf("LEG 2 ΔV: %.4f km s^-1\n", opt_dv_fb_norm);
fprintf("LEG 3 ΔV: %.4f km s^-1\n", opt_dv_rv_norm);

fprintf("\nOPTIMAL DATES\n")
fprintf("MERCURY DEP   MJD2000 %.3f\n", time_list(optimal_M_idx));
fprintf("MERCURY DEP   DATE    %.0f %.0f %.0f %.0f %.0f %.0f\n", mjd20002date(time_list(optimal_M_idx)));

fprintf("EARTH ARR/DEP MJD2000 %.3f\n", time_list(optimal_E_idx));
fprintf("EARTH ARR/DEP DATE    %.0f %.0f %.0f %.0f %.0f %.0f\n", mjd20002date(time_list(optimal_E_idx)));

fprintf("ASTEROID ARR  MJD2000 %.3f\n", time_list(optimal_A_idx));
fprintf("ASTEROID ARR  DATE    %.0f %.0f %.0f %.0f %.0f %.0f\n", mjd20002date(time_list(optimal_A_idx)));

% %% Plot the transfer trajectory for this mission
% % --- set up times and options ---
% mjd2k1 = time_list(optimal_M_idx);
% mjd2k2 = time_list(optimal_E_idx);
% mjd2k3 = time_list(optimal_A_idx);
% 
% t1 = mjd2k1 * 24 * 3600;
% t2 = mjd2k2 * 24 * 3600;
% t3 = mjd2k3 * 24 * 3600;
% t4 = t3 + (120*24*3600);
% t_full = linspace(t1, t4, 100);
% 
% options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);
% 
% % --- get planet and asteroid heliocentric positions and velocities ---
% [RM1, VM1] = get_planet_state(mjd2k1, planet_M_id, mu_sun); % mercury
% [RE1, VE1] = get_planet_state(mjd2k1, planet_E_id, mu_sun); % earth
% [RA1, VA1] = get_asteroid_state(mjd2k1, asteroid_id, mu_sun); % asteroid
% 
% % --- get satellite heliocentric positions and velocities ---
% RSM = RM1;
% VSM = reshape(V1_grid(optimal_M_idx, optimal_E_idx, :), [1, 3]); % at mercury
% 
% [RSE, ~] = get_planet_state(mjd2k2, 3, mu_sun);
% VSE = reshape(V3_grid(optimal_E_idx, optimal_A_idx, :), [1, 3]); % at earth (after GA)
% 
% [RSA, ~] = get_asteroid_state(mjd2k3, 316801, mu_sun);
% VSA = VA_list(optimal_A_idx, :); % at asteroid
% 
% % a. propagate the satellite's mercury-earth leg
% y_dep = [RSM, VSM];
% tspan_dep = t_full(t_full <= t2);
% [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_dep, y_dep, options);
% R_dep = Y(:, 1:3) ./ AU;
% 
% % b. propagate the satellite's earth-asteroid leg
% y_GA = [RSE, VSE];
% tspan_GA = t_full(t_full > t2 & t_full <= t3);
% [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_GA, y_GA, options);
% R_GA = Y(:, 1:3) ./ AU;
% 
% % c. propagate another month after asteroid rendez-vous
% y_rv = [RSA, VSA];
% tspan_rv = t_full(t_full > t3);
% [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_rv, y_rv, options);
% R_rv = Y(:, 1:3) ./ AU;
% 
% % d. propagate the departure planet orbit
% y_M = [RM1, VM1];
% [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), t_full, y_M, options);
% R_M = Y(:, 1:3) ./ AU;
% 
% % e. propagate the gravity-assist planet orbit
% y_E = [RE1, VE1];
% [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), t_full, y_E, options);
% R_E = Y(:, 1:3) ./ AU;
% 
% % f. propagate the arrival asteroid orbit
% y_A = [RA1, VA1];
% [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), t_full, y_A, options);
% R_A = Y(:, 1:3) ./ AU;
% 
% % --- plot ---
% figure("Name", "Orbit Plot"); hold on;
% % a. mercury-earth leg
% plot3(R_dep(:, 1), R_dep(:, 2), R_dep(:, 3), "y--");
% 
% % b. earth-asteroid leg
% plot3(R_GA(:, 1), R_GA(:, 2), R_GA(:, 3), "y");
% 
% % c. further asteroid propagation leg
% plot3(R_rv(:, 1), R_rv(:, 2), R_rv(:, 3), "y");
% 
% % d. planet mercury orbit
% plot3(R_M(:, 1), R_M(:, 2), R_M(:, 3), "r");
% 
% % e. planet earth orbit
% plot3(R_E(:, 1), R_E(:, 2), R_E(:, 3), "g");
% 
% % f. asteroid orbit
% plot3(R_A(:, 1), R_A(:, 2), R_A(:, 3), "b");
% 
% % g. planet mercury important positions (departure
% scatter3(R_M(1, 1), R_M(1, 2), R_M(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "r");
% 
% % h. planet earth important positions (departure and gravity assist)
% scatter3(R_E(1, 1), R_E(1, 2), R_E(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "g");
% scatter3(R_GA(1, 1), R_GA(1, 2), R_GA(1, 3), "filled", "MarkerFaceColor", "g");
% 
% % i. asteroid boundary positions (departure and rendez-vous)
% scatter3(R_A(1, 1), R_A(1, 2), R_A(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "b");
% scatter3(R_rv(1, 1), R_rv(1, 2), R_rv(1, 3), "filled", "MarkerFaceColor", "b");
% 
% % j. sun position
% scatter3(0, 0, 0, "filled", "MarkerFaceColor", "y");
% 
% % --- plot properties ---
% xlabel("X [AU]"); ylabel("Y [AU]"); zlabel("Z [AU]");
% title("Two-body problem orbit");
% legend( ...
%     "", ... a. mercury-earth leg
%     "", ... b. earth-asteroid leg
%     "", ... c. asteroid propagation
%     "", ... d. planet mercury orbit
%     "", ... e. planet earth orbit
%     "", ... f. asteroid orbit
%     "Mercury at departure", ...
%     "Earth at departure", ...
%     "Earth at GA", ...
%     "Asteroid at departure", ...
%     "Asteroid at rendez-vous", ...
%     "" ... j. sun position
%     );
% axis equal; grid on; view(3);
% hold off;
% 
% %% plot orbit animations
% % Define orbital data
% xM = R_M(:, 1); xE = R_E(:, 1); xA = R_A(:, 1); xS = [R_dep(:, 1); R_GA(:, 1); R_rv(:, 1);];
% yM = R_M(:, 2); yE = R_E(:, 2); yA = R_A(:, 2); yS = [R_dep(:, 2); R_GA(:, 2); R_rv(:, 2);];
% zM = R_M(:, 3); zE = R_E(:, 3); zA = R_A(:, 3); zS = [R_dep(:, 3); R_GA(:, 3); R_rv(:, 3);];
% 
% % Initialise the plot
% figure("Name", "Animated Orbit Plot"); hold on; grid on; view(3); axis equal;
% maxx = max([max(abs(xM)), max(abs(xE)), max(abs(xA)), max(abs(xS))]);
% maxy = max([max(abs(yM)), max(abs(yE)), max(abs(yA)), max(abs(yS))]);
% maxz = max([max(abs(zM)), max(abs(zE)), max(abs(zA)), max(abs(zS))]);
% xlim([-maxx, maxx]); ylim([-maxy, maxy]); zlim([-maxz, maxz])
% 
% scatter3(0, 0, 0, "filled", "MarkerFaceColor", "y");
% 
% % Create animated lines
% hM = animatedline("Color", "r", "MaximumNumPoints", inf);
% hE = animatedline("Color", "g", "MaximumNumPoints", inf);
% hA = animatedline("Color", "b", "MaximumNumPoints", inf);
% hS = animatedline("color", "y", "MaximumNumPoints", inf);
% 
% % Add markers for the "heads" of the comets
% headM = plot3(xM(1), yM(1), zM(1), "ro", "MarkerFaceColor", "r");
% headE = plot3(xE(1), yE(1), zE(1), "go", "MarkerFaceColor", "g");
% headA = plot3(xA(1), yA(1), zA(1), "bo", "MarkerFaceColor", "b");
% headS = plot3(xS(1), yS(1), zS(1), "yo", "MarkerFaceColor", "y");
% 
% % --- plot properties ---
% title("Simultaneous Multi-Orbit Animation");
% xlabel("X [AU]"); ylabel("Y [AU]"); zlabel("Z [AU]");
% legend( ...
%     "Sun Position", ...
%     "Mercury Orbit", ...
%     "Earth Orbit", ...
%     "Asteroid Orbit", ...
%     "Spacecraft Trajectory" ...
%     );
% 
% % Animation loop
% for i = 1:length(t_full)
%     % Update the tails
%     if i == 1
%         pause(2)
%     end
%     addpoints(hM, xM(i), yM(i), zM(i));
%     addpoints(hE, xE(i), yE(i), zE(i));
%     addpoints(hA, xA(i), yA(i), zA(i));
%     addpoints(hS, xS(i), yS(i), zS(i));
% 
%     % Update the heads
%     set(headM, "XData", xM(i), "YData", yM(i), "ZData", zM(i));
%     set(headE, "XData", xE(i), "YData", yE(i), "ZData", zE(i));
%     set(headA, "XData", xA(i), "YData", yA(i), "ZData", zA(i));
%     set(headS, "XData", xS(i), "YData", yS(i), "ZData", zS(i));
% 
%     drawnow limitrate; pause(0.03); 
% end
