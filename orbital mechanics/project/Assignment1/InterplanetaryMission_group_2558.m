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
% porkchop_plot( ...
%     "Mercury-Earth", ...
%     dv_grid1_norm, ...
%     time_list, ...
%     time_list);

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
% porkchop_plot( ...
%     "Earth-Asteroid", ...
%     dv_grid3_norm, ...
%     time_list, ...
%     time_list);

%% Manouvre 2: Gravity Assist at Earth
fprintf("\nconducting grid search 3 (flyby)... "); tic

possible_flyby_idxs = intersect(dv_grid1_valid_cols, dv_grid3_valid_rows);
rp_crit = planet_E_r + 500;

% Initialize global best tracker
opt_dv_tot = Inf;
opt_dv_launch = Inf; opt_dv_fb = Inf; opt_dv_rv = Inf;
opt_dv_launch_norm = Inf; opt_dv_fb_norm = Inf; opt_dv_rv_norm = Inf;
optimal_M_idx = NaN; optimal_E_idx = NaN; optimal_A_idx = NaN;

for j = 1:length(possible_flyby_idxs)
    jj = possible_flyby_idxs(j);
    
    % 1. Get Earth state and valid Mercury departures for this Earth date
    V_planet = VE_list(jj, :);
    valid_depart_idxs = find(~isnan(V2_grid(:, jj, 1))); % Indices of valid Mercury departures
    
    % 2. Pre-extract all valid Asteroid arrival velocities for this Earth date
    valid_arrive_idxs = dv_grid3_valid_cols; 
    
    V_plus_all = squeeze(V3_grid(jj, valid_arrive_idxs, :)); % all possible departures to asteroid
    
    % Pre-calculate v_inf_plus (departure from Earth)
    v_inf_plus_all = V_plus_all - V_planet;
    v_inf_plus_norms = vecnorm(v_inf_plus_all, 2, 2);
    
    % Pre-fetch Leg 3 costs
    dv_rv_costs = dv_grid3_norm(jj, valid_arrive_idxs)';
    
    % Loop over valid Mercury departures (i)
    for i = 1:length(valid_depart_idxs)
        ii = valid_depart_idxs(i);
        
        % --- HOISTED CALCULATIONS (Done once per i) ---
        V_minus = reshape(V2_grid(ii, jj, :), 1, 3);
        dv_launch_cost = dv_grid1_norm(ii, jj);
        
        % v_inf incoming (scalar for this i loop)
        v_inf_minus = V_minus - V_planet;
        v_inf_minus_norm = norm(v_inf_minus);
        
        % Calculate flyby limit criteria (scalar for this i loop)
        e_minus_crit = 1 + (rp_crit * v_inf_minus_norm^2) / planet_E_mu;
        delta_minus_term = 2 * asin(1 / e_minus_crit);
        
        % --- VECTORIZED OPERATIONS (Process all k at once) ---
        
        % 1. Calculate Flyby Delta V for all k
        % (V_plus_all is Nx3, V_minus is 1x3, automatic expansion handles this)
        dv_fb_vecs = V_plus_all - V_minus; 
        dv_fb_norms = vecnorm(dv_fb_vecs, 2, 2);
        
        % 2. Calculate Total Cost for all k
        dv_totals = dv_launch_cost + dv_fb_norms + dv_rv_costs;
        
        % 3. First Filter: Cost Limits (Boolean Mask)
        % Check if cost is valid AND better than global best
        mask_cost = (dv_fb_norms <= dv_lim) & (dv_totals < opt_dv_tot);
        
        if ~any(mask_cost)
            continue; % Skip expensive trig if no candidates passed
        end
        
        % --- EXPENSIVE CHECKS (Only on survivors) ---
        
        % Filter arrays to keep only candidates that passed cost check
        survivor_idxs = find(mask_cost);
        v_inf_plus_sub = v_inf_plus_all(survivor_idxs, :);
        v_inf_plus_norms_sub = v_inf_plus_norms(survivor_idxs);
        dv_totals_sub = dv_totals(survivor_idxs);
        dv_fb_norms_sub = dv_fb_norms(survivor_idxs);
        dv3_costs_sub = dv_rv_costs(survivor_idxs);
        
        % Calculate Turning Angle (Delta)
        % Dot product of single vector v_inf_minus with matrix v_inf_plus_sub
        dot_prods = v_inf_plus_sub * v_inf_minus';
        deltas = acos(dot_prods ./ (v_inf_plus_norms_sub * v_inf_minus_norm));
        
        % Calculate Max Turning Angle (Delta Crit)
        e_plus_crits = 1 + (rp_crit .* (v_inf_plus_norms_sub.^2)) ./ planet_E_mu;
        delta_plus_terms = 2 .* asin(1 ./ e_plus_crits);
        delta_crits = (delta_minus_term + delta_plus_terms) / 2;
        
        % 4. Final Filter: Geometry
        mask_geo = (deltas < delta_crits) & ~isnan(deltas);
        
        % --- UPDATE GLOBAL BEST ---
        if any(mask_geo)
            % Find the best among the survivors
            [min_val, best_idx_local] = min(dv_totals_sub(mask_geo));
            
            % If this batch contains a new global minimum
            if min_val < opt_dv_tot
                opt_dv_tot = min_val;
                
                % Map back to original indices
                valid_indices_subset = survivor_idxs(mask_geo);
                true_k_idx = valid_indices_subset(best_idx_local);
                
                optimal_M_idx = ii;
                optimal_E_idx = jj;
                optimal_A_idx = valid_arrive_idxs(true_k_idx);
                
                % Store specific components
                opt_dv_launch = dv_launch_cost;

                opt_dv_fb = dv_fb_norms_sub(mask_geo);
                opt_dv_fb = opt_dv_fb(best_idx_local);

                opt_dv_rv = dv3_costs_sub(mask_geo);
                opt_dv_rv = opt_dv_rv(best_idx_local);
            end
        end
    end
end
v_inf_minus_norm_best = norm(reshape(V2_grid(optimal_M_idx, optimal_E_idx, :), 1, 3)-VE_list(optimal_E_idx, :));
v_inf_plus_norm_best = norm(v_inf_plus_sub(best_idx_local));
eq = @(rp) deltas(best_idx_local) - ...
    asin(1 / (1 + (rp * v_inf_plus_norm_best^2) / planet_E_mu)) - ...
    asin(1 / (1 + (rp * v_inf_minus_norm_best^2) / planet_E_mu));
rp_ans = fzero(eq, rp_crit, optimset("Display", "off"));

disp("complete!");
toc

%% Final Results Output
fprintf("\nTOTAL ΔV REQUIRED: %.4f km s^-1\n", norm(opt_dv_tot));
fprintf("LEG 1 ΔV: %.4f km s^-1\n", norm(opt_dv_launch));
fprintf("LEG 2 ΔV: %.4f km s^-1\n", norm(opt_dv_fb));
fprintf("LEG 3 ΔV: %.4f km s^-1\n", norm(opt_dv_rv));

fprintf("\nOPTIMAL DATES\n")
fprintf("MERCURY DEP   MJD2000 %.3f\n", time_list(optimal_M_idx));
fprintf("MERCURY DEP   DATE    %.0f %.0f %.0f %.0f %.0f %.0f\n", mjd20002date(time_list(optimal_M_idx)));

fprintf("EARTH ARR/DEP MJD2000 %.3f\n", time_list(optimal_E_idx));
fprintf("EARTH ARR/DEP DATE    %.0f %.0f %.0f %.0f %.0f %.0f\n", mjd20002date(time_list(optimal_E_idx)));

fprintf("ASTEROID ARR  MJD2000 %.3f\n", time_list(optimal_A_idx));
fprintf("ASTEROID ARR  DATE    %.0f %.0f %.0f %.0f %.0f %.0f\n", mjd20002date(time_list(optimal_A_idx)));

%% Manouvre 2: Gravity Assist at Earth (Fully Vectorized)
fprintf("\nconducting grid search 3 (flyby)... "); tic

possible_flyby_idxs = intersect(dv_grid1_valid_cols, dv_grid3_valid_rows);
rp_crit = planet_E_r + 500;

% Global best tracker
lowest_dvtot_an = Inf;
opt_dv1 = NaN; opt_dv2 = NaN; opt_dv3 = NaN;
optimal_M_idx = NaN; optimal_E_idx = NaN; optimal_A_idx = NaN;

% Pre-extract Earth Velocity matrix
VE_matrix = VE_list; 

for j = 1:length(possible_flyby_idxs)
    jj = possible_flyby_idxs(j);
    
    % --- 1. PREPARE DATA VECTORS ---
    V_planet = VE_matrix(jj, :);
    
    % A. Mercury Departures (Rows)
    valid_depart_idxs = find(~isnan(V2_grid(:, jj, 1))); 
    if isempty(valid_depart_idxs), continue; end
    
    % Extract vectors (N_dep x 3)
    V_minus_all = squeeze(V2_grid(valid_depart_idxs, jj, :));
    if size(V_minus_all, 2) == 1, V_minus_all = V_minus_all'; end 
    
    % Extract costs (N_dep x 1)
    dv1_costs = dv_grid1_norm(valid_depart_idxs, jj);

    % B. Asteroid Arrivals (Columns)
    valid_arrive_idxs = dv_grid3_valid_cols; 
    
    % Extract vectors (N_arr x 3)
    V_plus_all = squeeze(V3_grid(jj, valid_arrive_idxs, :));
    if size(V_plus_all, 2) == 1, V_plus_all = V_plus_all'; end
    
    % Extract costs (1 x N_arr) - Note the transpose to make it a row vector
    dv3_costs = dv_grid3_norm(jj, valid_arrive_idxs)';
    
    % --- 2. VECTORIZED CALCULATIONS (Implicit Expansion) ---
    
    % Reshape for broadcasting:
    % V_minus: (N_dep x 1 x 3)
    % V_plus:  (1 x N_arr x 3)
    V_minus_b = permute(V_minus_all, [1 3 2]); 
    V_plus_b  = permute(V_plus_all, [3 1 2]);  
    
    % Calculate Flyby Delta V Matrix (N_dep x N_arr)
    % MATLAB automatically expands dimensions to subtract every row from every col
    dv_fb_vecs = V_plus_b - V_minus_b; 
    dv_fb_mat = vecnorm(dv_fb_vecs, 2, 3);
    
    % Calculate Total Cost Matrix (N_dep x N_arr)
    % (N_dep x 1) + (N_dep x N_arr) + (1 x N_arr)
    dv_total_mat = dv1_costs + dv_fb_mat + dv3_costs;
    
    % --- 3. FILTERING (Boolean Masks) ---
    
    % Filter 1: Cost thresholds
    mask_cost = (dv_fb_mat <= dv_lim) & (dv_total_mat < lowest_dvtot_an);
    
    if ~any(mask_cost, 'all')
        continue; 
    end
    
    % --- 4. GEOMETRY CHECKS (Vectorized) ---
    
    % We calculate geometry ONLY for candidates passing the cost filter.
    % To do this efficiently, we use linear indexing on the survivors.
    idx_survivors = find(mask_cost);
    
    % Map linear indices back to Row (Mercury) and Col (Asteroid) subscripts
    [sub_row, sub_col] = ind2sub(size(mask_cost), idx_survivors);
    
    % Extract relevant velocity vectors for survivors
    V_min_surv = V_minus_all(sub_row, :);
    V_plus_surv = V_plus_all(sub_col, :);
    
    % Calculate v_inf vectors
    v_inf_minus = V_min_surv - V_planet;
    v_inf_plus  = V_plus_surv - V_planet;
    
    vn_minus = vecnorm(v_inf_minus, 2, 2);
    vn_plus  = vecnorm(v_inf_plus, 2, 2);
    
    % Turning Angle (Delta)
    dot_prod = dot(v_inf_minus, v_inf_plus, 2);
    deltas = acos(dot_prod ./ (vn_minus .* vn_plus));
    
    % Max Turning Angle (Delta Crit)
    e_minus = 1 + (rp_crit .* vn_minus.^2) ./ planet_E_mu;
    e_plus  = 1 + (rp_crit .* vn_plus.^2) ./ planet_E_mu;
    
    delta_crit = asin(1./e_minus) + asin(1./e_plus); % (2*asin...)/2 simplifies
    
    % Filter 2: Geometry
    mask_geo = (deltas < delta_crit) & ~isnan(deltas);
    
    if any(mask_geo)
        % Get costs of survivors that passed geometry
        survivor_costs = dv_total_mat(idx_survivors(mask_geo));
        
        [min_val, best_loc_idx] = min(survivor_costs);
        
        if min_val < lowest_dvtot_an
            lowest_dvtot_an = min_val;
            
            % Trace back the indices
            % 1. Which survivor was it?
            best_linear_idx_in_grid = idx_survivors(mask_geo);
            best_linear_idx_in_grid = best_linear_idx_in_grid(best_loc_idx);
            
            % 2. Convert grid index to row/col
            [r, c] = ind2sub(size(dv_total_mat), best_linear_idx_in_grid);
            
            % 3. Map to real global indices
            optimal_M_idx = valid_depart_idxs(r);
            optimal_E_idx = jj;
            optimal_A_idx = valid_arrive_idxs(c);
            
            % Store components
            opt_dv1 = dv1_costs(r);
            opt_dv2 = dv_fb_mat(best_linear_idx_in_grid);
            opt_dv3 = dv3_costs(c);
        end
    end
end
disp("complete!");
toc

%% Final Results Output
fprintf("\nTOTAL ΔV REQUIRED: %.4f km s^-1\n", norm(lowest_dvtot_an));
fprintf("LEG 1 ΔV: %.4f km s^-1\n", norm(opt_dv1));
fprintf("LEG 2 ΔV: %.4f km s^-1\n", norm(opt_dv2));
fprintf("LEG 3 ΔV: %.4f km s^-1\n", norm(opt_dv3));

fprintf("\nOPTIMAL DATES\n")
fprintf("MERCURY DEP   MJD2000 %.3f\n", time_list(optimal_M_idx));
fprintf("MERCURY DEP   DATE    %.0f %.0f %.0f %.0f %.0f %.0f\n", mjd20002date(time_list(optimal_M_idx)));

fprintf("EARTH ARR/DEP MJD2000 %.3f\n", time_list(optimal_E_idx));
fprintf("EARTH ARR/DEP DATE    %.0f %.0f %.0f %.0f %.0f %.0f\n", mjd20002date(time_list(optimal_E_idx)));

fprintf("ASTEROID ARR  MJD2000 %.3f\n", time_list(optimal_A_idx));
fprintf("ASTEROID ARR  DATE    %.0f %.0f %.0f %.0f %.0f %.0f\n", mjd20002date(time_list(optimal_A_idx)));

%% Plot the transfer trajectory for this mission
% --- set up times and options ---
mjd2k1 = time_list(optimal_M_idx);
mjd2k2 = time_list(optimal_E_idx);
mjd2k3 = time_list(optimal_A_idx);

t1 = mjd2k1 * 24 * 3600;
t2 = mjd2k2 * 24 * 3600;
t3 = mjd2k3 * 24 * 3600;
t4 = t3 + (120*24*3600);
t_full = linspace(t1, t4, 500);

options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% --- get planet and asteroid heliocentric positions and velocities ---
[RM1, VM1] = get_planet_state(mjd2k1, planet_M_id, mu_sun); % mercury
[RE1, VE1] = get_planet_state(mjd2k1, planet_E_id, mu_sun); % earth
[RA1, VA1] = get_asteroid_state(mjd2k1, asteroid_id, mu_sun); % asteroid

% --- get satellite heliocentric positions and velocities ---
RSM = RM1;
VSM = reshape(V1_grid(optimal_M_idx, optimal_E_idx, :), [1, 3]); % at mercury

[RSE, ~] = get_planet_state(mjd2k2, 3, mu_sun);
VSE = reshape(V3_grid(optimal_E_idx, optimal_A_idx, :), [1, 3]); % at earth (after GA)

[RSA, ~] = get_asteroid_state(mjd2k3, 316801, mu_sun);
VSA = VA_list(optimal_A_idx, :); % at asteroid

% a. propagate the satellite's mercury-earth leg
y_dep = [RSM, VSM];
tspan_dep = t_full(t_full <= t2);
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_dep, y_dep, options);
R_dep = Y(:, 1:3) ./ AU;

% b. propagate the satellite's earth-asteroid leg
y_GA = [RSE, VSE];
tspan_GA = t_full(t_full > t2 & t_full <= t3);
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_GA, y_GA, options);
R_GA = Y(:, 1:3) ./ AU;

% c. propagate another month after asteroid rendez-vous
y_rv = [RSA, VSA];
tspan_rv = t_full(t_full > t3);
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), tspan_rv, y_rv, options);
R_rv = Y(:, 1:3) ./ AU;

% d. propagate the departure planet orbit
y_M = [RM1, VM1];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), t_full, y_M, options);
R_M = Y(:, 1:3) ./ AU;

% e. propagate the gravity-assist planet orbit
y_E = [RE1, VE1];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), t_full, y_E, options);
R_E = Y(:, 1:3) ./ AU;

% f. propagate the arrival asteroid orbit
y_A = [RA1, VA1];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), t_full, y_A, options);
R_A = Y(:, 1:3) ./ AU;

% --- plot ---
figure("Name", "Orbit Plot"); hold on;
% a. mercury-earth leg
plot3(R_dep(:, 1), R_dep(:, 2), R_dep(:, 3), "y--");

% b. earth-asteroid leg
plot3(R_GA(:, 1), R_GA(:, 2), R_GA(:, 3), "y");

% c. further asteroid propagation leg
plot3(R_rv(:, 1), R_rv(:, 2), R_rv(:, 3), "y");

% d. planet mercury orbit
plot3(R_M(:, 1), R_M(:, 2), R_M(:, 3), "r");

% e. planet earth orbit
plot3(R_E(:, 1), R_E(:, 2), R_E(:, 3), "g");

% f. asteroid orbit
plot3(R_A(:, 1), R_A(:, 2), R_A(:, 3), "b");

% g. planet mercury important positions (departure
scatter3(R_M(1, 1), R_M(1, 2), R_M(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "r");

% h. planet earth important positions (departure and gravity assist)
scatter3(R_E(1, 1), R_E(1, 2), R_E(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "g");
scatter3(R_GA(1, 1), R_GA(1, 2), R_GA(1, 3), "filled", "MarkerFaceColor", "g");

% i. asteroid boundary positions (departure and rendez-vous)
scatter3(R_A(1, 1), R_A(1, 2), R_A(1, 3), "MarkerFaceColor", "none", "MarkerEdgeColor", "b");
scatter3(R_rv(1, 1), R_rv(1, 2), R_rv(1, 3), "filled", "MarkerFaceColor", "b");

% j. sun position
scatter3(0, 0, 0, "filled", "MarkerFaceColor", "y");

% --- plot properties ---
xlabel("X [AU]"); ylabel("Y [AU]"); zlabel("Z [AU]");
title("Two-body problem orbit");
legend( ...
    "", ... a. mercury-earth leg
    "", ... b. earth-asteroid leg
    "", ... c. asteroid propagation
    "", ... d. planet mercury orbit
    "", ... e. planet earth orbit
    "", ... f. asteroid orbit
    "Mercury at departure", ...
    "Earth at departure", ...
    "Earth at GA", ...
    "Asteroid at departure", ...
    "Asteroid at rendez-vous", ...
    "" ... j. sun position
    );
axis equal; grid on; view(3);
hold off;

%% plot orbit animations
% Define orbital data
xM = R_M(:, 1); xE = R_E(:, 1); xA = R_A(:, 1);
yM = R_M(:, 2); yE = R_E(:, 2); yA = R_A(:, 2);
zM = R_M(:, 3); zE = R_E(:, 3); zA = R_A(:, 3);

xS = [R_dep(:, 1); R_GA(:, 1); R_rv(:, 1);];
yS = [R_dep(:, 2); R_GA(:, 2); R_rv(:, 2);];
zS = [R_dep(:, 3); R_GA(:, 3); R_rv(:, 3);];

% Initialise the plot
figure("Name", "Animated Orbit Plot"); hold on; grid on; view(3); axis equal;
maxx = max([max(abs(xM)), max(abs(xE)), max(abs(xA)), max(abs(xS))]);
maxy = max([max(abs(yM)), max(abs(yE)), max(abs(yA)), max(abs(yS))]);
maxz = max([max(abs(zM)), max(abs(zE)), max(abs(zA)), max(abs(zS))]);
xlim([-maxx, maxx]); ylim([-maxy, maxy]); zlim([-maxz, maxz])

scatter3(0, 0, 0, "filled", "MarkerFaceColor", "y");

% Create animated lines
hM = animatedline("Color", "r", "MaximumNumPoints", inf);
hE = animatedline("Color", "g", "MaximumNumPoints", inf);
hA = animatedline("Color", "b", "MaximumNumPoints", inf);
hS = animatedline("color", "y", "MaximumNumPoints", inf);

% Add markers for the "heads" of the comets
headM = plot3(xM(1), yM(1), zM(1), "ro", "MarkerFaceColor", "r");
headE = plot3(xE(1), yE(1), zE(1), "go", "MarkerFaceColor", "g");
headA = plot3(xA(1), yA(1), zA(1), "bo", "MarkerFaceColor", "b");
headS = plot3(xS(1), yS(1), zS(1), "yo", "MarkerFaceColor", "y");

% --- plot properties ---
title("Simultaneous Multi-Orbit Animation");
xlabel("X [AU]"); ylabel("Y [AU]"); zlabel("Z [AU]");
legend( ...
    "Sun Position", ...
    "Mercury Orbit", ...
    "Earth Orbit", ...
    "Asteroid Orbit", ...
    "Spacecraft Trajectory" ...
    );

% Animation loop
for i = 1:length(t_full)
    % Update the tails
    if i == 1
        pause(2)
    end
    addpoints(hM, xM(i), yM(i), zM(i));
    addpoints(hE, xE(i), yE(i), zE(i));
    addpoints(hA, xA(i), yA(i), zA(i));
    addpoints(hS, xS(i), yS(i), zS(i));

    % Update the heads
    set(headM, "XData", xM(i), "YData", yM(i), "ZData", zM(i));
    set(headE, "XData", xE(i), "YData", yE(i), "ZData", zE(i));
    set(headA, "XData", xA(i), "YData", yA(i), "ZData", zA(i));
    set(headS, "XData", xS(i), "YData", yS(i), "ZData", zS(i));

    drawnow limitrate; pause(0.03); 
end
