%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% ASSIGNMENT 1: INTERPLANETARY EXPLORER MISSION %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Departure planet: Mercury
% Flyby planet: Earth
% Arrival body: Asteroid N.316801
% Earliest departure: 00:00:00 01/01/2030
% Latest arrival: 00:00:00 01/01/2060

clear
close all
clc

%% Data
mu_sun = astroConstants(4);
AU = astroConstants(2);
mu_earth = astroConstants(13);
R_earth = astroConstants(23);
h_atm = 300;                           % height of the Earth's atmosphere
safety = 200;                          % safety margin
r_p_min = R_earth + h_atm + safety;    % minimum radius for the fly-by

%% Estimates of the first heliocentric leg: synodic period, time of flight
% Important assumptions: orbits are coplanar and circular

% Average distance of Mercury and Earth from the Sun
r_M = 0.3871 * AU;
r_E = 1.0000 * AU;

% Orbital period
T_M = 2 * pi * sqrt(r_M^3 / mu_sun);
T_E = 2 * pi * sqrt(r_E^3 / mu_sun);

% Synodic period
T_syn1 = T_M * T_E / abs(T_M - T_E);

% Time of flight
a1_H = (r_M + r_E) / 2;
T1_H = 2 * pi * sqrt(a1_H^3 / mu_sun);
TOF1 = T1_H / 2;

fprintf('---- ESTIMATES FIRST HELIOCENTRIC LEG (MERCURY-EARTH) ---- \n');
fprintf('Synodic period = %.2f [days]\n', T_syn1 / 86400);
fprintf('TOF = %.2f [days] \n\n', TOF1 / 86400);

%% Estimates of the second heliocentric leg: synodic period, time of flight
date_window = [2030, 1, 1, 0, 0, 0;
               2060, 1, 1, 0, 0, 0];
dates = dateWindow(date_window, 1);

N = size(dates, 1);
dates_mjd2000 = zeros(N, 1);
kepAsteroid = zeros(N, 6);

for ii = 1:N
    dates_mjd2000(ii) = date2mjd2000(dates(ii, :));
    [kep, ~, ~] = ephAsteroids(dates_mjd2000(ii), 316801);
    kepAsteroid(ii,:) = kep;
end

% Averaged values
aAst_averaged = mean(kepAsteroid(:,1));
eAst_averaged = mean(kepAsteroid(:,2));
iAst_averaged = mean(kepAsteroid(:,3));
OMAst_averaged = mean(kepAsteroid(:,4));
omAst_averaged = mean(kepAsteroid(:,5));

% Average distance of Asteroid 316801 from the Sun
r_A = aAst_averaged;
T_A = 2 * pi * sqrt(r_A^3 / mu_sun);

% Synodic period
T_syn2 = T_E * T_A / abs(T_E - T_A);

% Time of flight
a2_H = (r_E + r_A) / 2;
T2_H = 2*pi*sqrt(a2_H^3 / mu_sun);
TOF2 = T2_H / 2;

fprintf('---- ESTIMATES SECOND HELIOCENTRIC LEG (EARTH-ASTEROID) ---- \n');
fprintf('Synodic period = %.2f [days]\n', T_syn2 / 86400);
fprintf('TOF = %.2f [days] \n\n', TOF2 / 86400);

%% Values for subsequent calculations based on preliminary estimates
% First heliocentric leg
T_syn_ME = 116;    % days
TOF_ME = 105;      % days
TOF_ME_min = 80;   % days
TOF_ME_max = 130;  % days

% Second heliocentric leg
T_syn_EA = 472;    % days
TOF_EA = 459;      % days
TOF_EA_min = 300;  % days
TOF_EA_max = 650;  % days

%% Analysis and porkchop plot for the first synodic period (first leg)
t0 = date2mjd2000([2030 1 1 0 0 0]);
dep_dates_ME = t0 : 1 : t0 + T_syn_ME;
arr_dates_ME = dep_dates_ME(1) + TOF_ME_min : 1 : dep_dates_ME(end) + TOF_ME_max;

DV_matrix_ME = zeros(length(dep_dates_ME), length(arr_dates_ME));
TOF_matrix_ME = zeros(length(dep_dates_ME), length(arr_dates_ME));
DV_min_ME = 1e4;

for ii = 1:length(dep_dates_ME)
    for jj = 1:length(arr_dates_ME)
        if arr_dates_ME(jj) > dep_dates_ME(ii)
            
            % Check for valid TOF
            TOF_check = arr_dates_ME(jj) - dep_dates_ME(ii);
            if TOF_check < TOF_ME_min || TOF_check > TOF_ME_max
                DV_matrix_ME(ii,jj) = NaN;
                continue
            end
            
            [kep_dep, ~] = uplanet(dep_dates_ME(ii), 1);
            [r1, v1] = par2car(kep_dep(1), kep_dep(2), kep_dep(3), ...
                              kep_dep(4), kep_dep(5), kep_dep(6), mu_sun);
            
            [kep_arr, ~] = uplanet(arr_dates_ME(jj), 3);
            [r2, v2] = par2car(kep_arr(1), kep_arr(2), kep_arr(3), ...
                              kep_arr(4), kep_arr(5), kep_arr(6), mu_sun);
            
            TOF_ME = (arr_dates_ME(jj) - dep_dates_ME(ii)) * 86400;
            [~, ~, ~, ~, v_dep, v_arr, ~, ~] = lambertMR(r1, r2, TOF_ME, mu_sun);
            
            DV_matrix_ME(ii,jj) = norm(v_dep' - v1) + norm(v_arr' - v2);
            TOF_matrix_ME(ii,jj) = TOF_ME / 86400;
        else
            DV_matrix_ME(ii,jj) = NaN;
        end

        % Store the best results
        if DV_matrix_ME(ii,jj) < DV_min_ME
            DV_min_ME = DV_matrix_ME(ii,jj);
            optimal_TOF_ME = TOF_matrix_ME(ii,jj);
            optimal_dep_date_ME = dep_dates_ME(ii);
            optimal_arr_date_ME = arr_dates_ME(jj);
        end
    end
end

% Optimal departure and arrival dates
optimal_dep_date_datenum_ME = datenum(mjd20002date(optimal_dep_date_ME));
optimal_arr_date_datenum_ME = datenum(mjd20002date(optimal_arr_date_ME));
optimal_dep_date_string_ME = datestr(optimal_dep_date_datenum_ME, 'dd/mm/yyyy');
optimal_arr_date_string_ME = datestr(optimal_arr_date_datenum_ME, 'dd/mm/yyyy');

fprintf('---- DELTA_V MINIMUM (MERCURY-EARTH) ---- \n');
fprintf('DV_min = %.4f [km/s] \n', DV_min_ME);
fprintf('Optimal TOF = %d [days] \n', optimal_TOF_ME);
fprintf('Optimal departure date: %s \n', optimal_dep_date_string_ME);
fprintf('Optimal arrival date: %s \n\n', optimal_arr_date_string_ME);

% Porkchop plot
dep_datenum_ME = zeros(size(dep_dates_ME));
for kk = 1:length(dep_dates_ME)
    dep_datenum_ME(kk) = datenum(mjd20002date(dep_dates_ME(kk)));
end

arr_datenum_ME = zeros(size(arr_dates_ME));
for kk = 1:length(arr_dates_ME)
    arr_datenum_ME(kk) = datenum(mjd20002date(arr_dates_ME(kk)));
end

[departure_grid_ME, arrival_grid_ME] = meshgrid(dep_datenum_ME, arr_datenum_ME);
DV_plot_ME = DV_matrix_ME';

figure()
[C,h] = contour(departure_grid_ME, arrival_grid_ME, DV_plot_ME, 'LineWidth', 1.2);
clabel(C,h,'labelspacing',300);
clim([20 120]);
cb = colorbar;
colormap('parula')
cmap = colormap;
n_colors = size(cmap, 1);

color_index = round((DV_min_ME - 20) / (120 - 20) * (n_colors - 1)) + 1;
color_index = max(1, min(n_colors, color_index));
optimal_color = cmap(color_index, :);

hold on;
plot(optimal_dep_date_datenum_ME, optimal_arr_date_datenum_ME, 'o', ...
'MarkerFaceColor', optimal_color, 'MarkerEdgeColor', 'k', 'MarkerSize', 6);

% Selection (and plot) of the windows limits of interest
dep_minus_20 = optimal_dep_date_datenum_ME - 20;
dep_plus_20 = optimal_dep_date_datenum_ME + 20;
y_limits = ylim;
plot([dep_minus_20 dep_minus_20], y_limits, '--k', 'LineWidth', 1.5);
plot([dep_plus_20 dep_plus_20], y_limits, '--k', 'LineWidth', 1.5);
datetick('x','yyyy-mmm-dd','keeplimits')
datetick('y','yyyy-mmm-dd','keeplimits')
x_limits = xlim;
y_limits = ylim;
set(gca, 'XTick', linspace(x_limits(1), x_limits(2), 10));
set(gca, 'YTick', linspace(y_limits(1), y_limits(2), 10));

datetick('x','yyyy-mmm-dd','keeplimits','keepticks')
datetick('y','yyyy-mmm-dd','keeplimits','keepticks')
xlabel('Departure date')
ylabel('Arrival date')
ylabel(cb, '\DeltaV [km/s]')
title('Porkchop plot Mercury - Earth')
grid on;

%% Analysis and porkchop plot for the first synodic period (second leg)
t0 = date2mjd2000([2030 1 1 0 0 0]);
dep_dates_EA = t0 : 1 : t0 + T_syn_EA;
arr_dates_EA = dep_dates_EA(1) + TOF_EA_min : 5 : dep_dates_EA(end) + TOF_EA_max;

DV_matrix_EA = zeros(length(dep_dates_EA), length(arr_dates_EA));
TOF_matrix_EA = zeros(length(dep_dates_EA), length(arr_dates_EA));
DV_min_EA = 1e4;

for ii = 1:length(dep_dates_EA)
    for jj = 1:length(arr_dates_EA)
        if arr_dates_EA(jj) > dep_dates_EA(ii)
            
            % Check for valid TOF
            TOF_check = arr_dates_EA(jj) - dep_dates_EA(ii);
            if TOF_check < TOF_EA_min || TOF_check > TOF_EA_max
                DV_matrix_EA(ii,jj) = NaN;
                continue
            end
            
            [kep_dep, ~] = uplanet(dep_dates_EA(ii), 3);
            [r1, v1] = par2car(kep_dep(1), kep_dep(2), kep_dep(3), ...
                              kep_dep(4), kep_dep(5), kep_dep(6), mu_sun);
            
            [kep_arr, ~, ~] = ephAsteroids(arr_dates_EA(jj), 316801);
            [r2, v2] = par2car(kep_arr(1), kep_arr(2), kep_arr(3), ...
                              kep_arr(4), kep_arr(5), kep_arr(6), mu_sun);
            
            TOF_EA = (arr_dates_EA(jj) - dep_dates_EA(ii)) * 86400;
            [~, ~, ~, ~, v_dep, v_arr, ~, ~] = lambertMR(r1, r2, TOF_EA, mu_sun);
            
            DV_matrix_EA(ii,jj) = norm(v_dep' - v1) + norm(v_arr' - v2);
            TOF_matrix_EA(ii,jj) = TOF_EA / 86400;
        else
            DV_matrix_EA(ii,jj) = NaN;
        end

        % Store the best results
        if DV_matrix_EA(ii,jj) < DV_min_EA
            DV_min_EA = DV_matrix_EA(ii,jj);
            optimal_TOF_EA = TOF_matrix_EA(ii,jj);
            optimal_dep_date_EA = dep_dates_EA(ii);
            optimal_arr_date_EA = arr_dates_EA(jj);
        end
    end
end

% Optimal departure and arrival dates
optimal_dep_date_datenum_EA = datenum(mjd20002date(optimal_dep_date_EA));
optimal_arr_date_datenum_EA = datenum(mjd20002date(optimal_arr_date_EA));
optimal_dep_date_string_EA = datestr(optimal_dep_date_datenum_EA, 'dd/mm/yyyy');
optimal_arr_date_string_EA = datestr(optimal_arr_date_datenum_EA, 'dd/mm/yyyy');

fprintf('---- DELTA_V MINIMUM (EARTH-ASTEROID) ---- \n');
fprintf('DV_min = %.4f [km/s] \n', DV_min_EA);
fprintf('Optimal TOF = %d [days] \n', optimal_TOF_EA);
fprintf('Optimal departure date: %s \n', optimal_dep_date_string_EA);
fprintf('Optimal arrival date: %s \n\n', optimal_arr_date_string_EA);

% Porkchop plot
dep_datenum_EA = zeros(size(dep_dates_EA));
for kk = 1:length(dep_dates_EA)
    dep_datenum_EA(kk) = datenum(mjd20002date(dep_dates_EA(kk)));
end

arr_datenum_EA = zeros(size(arr_dates_EA));
for kk = 1:length(arr_dates_EA)
    arr_datenum_EA(kk) = datenum(mjd20002date(arr_dates_EA(kk)));
end

[departure_grid_EA, arrival_grid_EA] = meshgrid(dep_datenum_EA, arr_datenum_EA);
DV_plot_EA = DV_matrix_EA';

figure()
[C,h] = contour(departure_grid_EA, arrival_grid_EA, DV_plot_EA, 'LineWidth', 1.2);
clabel(C,h,'labelspacing',300);
clim([10 60]);
cb = colorbar;
colormap('parula')
cmap = colormap;
n_colors = size(cmap, 1);

color_index = round((DV_min_EA - 10) / (60 - 10) * (n_colors - 1)) + 1;
color_index = max(1, min(n_colors, color_index));
optimal_color = cmap(color_index, :);

hold on;
plot(optimal_dep_date_datenum_EA, optimal_arr_date_datenum_EA, 'o', ...
'MarkerFaceColor', optimal_color, 'MarkerEdgeColor', 'k', 'MarkerSize', 6);

% Selection (and plot) of the windows limits of interest
dep_minus_40 = optimal_dep_date_datenum_EA - 40;
dep_plus_40 = optimal_dep_date_datenum_EA + 40;
y_limits = ylim;
plot([dep_minus_40 dep_minus_40], y_limits, '--k', 'LineWidth', 1.5);
plot([dep_plus_40 dep_plus_40], y_limits, '--k', 'LineWidth', 1.5);
datetick('x','yyyy-mmm-dd','keeplimits')
datetick('y','yyyy-mmm-dd','keeplimits')
x_limits = xlim;
y_limits = ylim;
set(gca, 'XTick', linspace(x_limits(1), x_limits(2), 10));
set(gca, 'YTick', linspace(y_limits(1), y_limits(2), 10));

datetick('x','yyyy-mmm-dd','keeplimits','keepticks')
datetick('y','yyyy-mmm-dd','keeplimits','keepticks')
xlabel('Departure date')
ylabel('Arrival date')
ylabel(cb, '\DeltaV [km/s]')
title('Porkchop plot Earth - Asteroid')
grid on;

%% Windows of interest and intersections
% Latest possible arrival date at the asteroid
t_max = date2mjd2000([2060 1 1 0 0 0]);

% First departure window of interest from Mercury
t_dep1_0 = date2mjd2000([2030 4 6 0 0 0]);
t_dep1 = -20 + t_dep1_0 : 1 : t_dep1_0 + 20;

T_syn1 = 116;
TOF1_min = 80;
TOF1_max = 130;

% First departure window of interest from Earth
t_dep2_0 = date2mjd2000([2030 7 12 0 0 0]);
t_dep2 = -40 + t_dep2_0 : 1 : t_dep2_0 + 40;
T_syn2 = 472;
TOF2_min = 300;
TOF2_max = 650;

% Avaiable windows of interest between 01/01/2030 and 01/01/2060
% worst scenario: t_dep1(end) + K1*T_syn_ME + TOF_ME_max + TOF_EA_max <= t_max
K1 = floor((t_max - t_dep1(end) - TOF1_max - TOF2_max) / T_syn1);

% worst scenario: t_dep2(end) + K2*T_syn_EA + TOF_EA_max <= t_max
K2 = floor((t_max - t_dep2(end) - TOF2_max) / T_syn2);

fprintf('---- VALUES FOR K1 AND K2 ---- \n');
fprintf('K1 = %d (number of available windows: %d) \n', K1, K1 + 1);
fprintf('K2 = %d (number of available windows: %d) \n\n', K2, K2 + 1);

% All possible departure dates from Mercury
t_dep1_all = [];
for k = 0:K1
    t_dep1_k = t_dep1 + k*T_syn1;
    t_dep1_all = [t_dep1_all; t_dep1_k];
end

% All possible departure dates from Earth
t_dep2_all = [];
for k = 0:K2
    t_dep2_k = t_dep2 + k*T_syn2;
    t_dep2_all = [t_dep2_all; t_dep2_k];
end

% Vectors for possible TOFs
TOF1_vec = TOF1_min : 1 : TOF1_max;
TOF2_vec = TOF2_min : 2 : TOF2_max;

% All possible arrival dates at Earth
t_arr1_all = zeros(size(t_dep1_all, 1), size(t_dep1_all, 2), length(TOF1_vec));

for ii = 1:size(t_dep1_all, 1)
    for jj = 1:size(t_dep1_all, 2)
        for kk = 1:length(TOF1_vec)
            t_arr1_all(ii, jj, kk) = t_dep1_all(ii, jj) + TOF1_vec(kk);
        end
    end
end

% All possible arrival dates at Asteroid 316801
t_arr2_all = zeros(size(t_dep2_all, 1), size(t_dep2_all, 2), length(TOF2_vec));

for ii = 1:size(t_dep2_all, 1)
    for jj = 1:size(t_dep2_all, 2)
        for kk = 1:length(TOF2_vec)
            t_arr2_all(ii, jj, kk) = t_dep2_all(ii, jj) + TOF2_vec(kk);
        end
    end
end

% Intersection to find all the possible flyby dates
t_arr1_vec = t_arr1_all(:);
t_dep2_vec = t_dep2_all(:);
[t_GA_common, idx_arr1, idx_dep2] = intersect(t_arr1_vec, t_dep2_vec);

fprintf('---- FLYBY DATES ---- \n');
fprintf('Number of exact flyby dates: %d \n', length(t_GA_common));

% Understanding which window each possible Earth flyby date belongs to
n_windows_1 = size(t_arr1_all, 1); % windows Mercury - Earth
n_days_1 = size(t_arr1_all, 2); % days in each window
n_tof_1 = size(t_arr1_all, 3); % TOF for each day in each window

n_windows_2 = size(t_dep2_all, 1); % windows Earth - Asteroid
n_days_2 = size(t_dep2_all, 2); % days in each window
n_tof_2 = size(t_dep2_all, 3); % TOF for each day in each window

% All possible flyby dates
K1_indices = zeros(length(t_GA_common), 1);
K2_indices = zeros(length(t_GA_common), 1);

for ii = 1:length(t_GA_common)

    % conversion from linear index to matrix index
    % k1 is the index for the synodic window Mercury - Earth
    [k1, ~, ~] = ind2sub([n_windows_1, n_days_1, n_tof_1], idx_arr1(ii));
    K1_indices(ii) = k1 - 1;
    
    % k2 is the index for the synodic window Earth - Asteroid
    [k2, ~, ~] = ind2sub([n_windows_2, n_days_2, n_tof_2], idx_dep2(ii));
    K2_indices(ii) = k2 - 1;
end

% Mask to understand if the flyby is valid or to be discarded
mask = true(length(t_GA_common), 1);

for ii = 1:length(t_GA_common)
    t_GA = t_GA_common(ii);
    
    % Quick estimate with average TOF of the window, to see if it is worth
    % doing the more expensive grid search
    t_dep1_est = t_GA - (TOF1_min + TOF1_max)/2;
    t_arr2_est = t_GA + (TOF2_min + TOF2_max)/2;
    
    [kepM, ~] = uplanet(t_dep1_est, 1);
    [rM, vM] = par2car(kepM(1), kepM(2), kepM(3), ...
                           kepM(4), kepM(5), kepM(6), mu_sun);

    [kepE, ~] = uplanet(t_GA, 3);
    [rE, vE] = par2car(kepE(1), kepE(2), kepE(3), ...
                           kepE(4), kepE(5), kepE(6), mu_sun);

    [kepA,~,~] = ephAsteroids(t_arr2_est, 316801);
    [rA, vA] = par2car(kepA(1), kepA(2), kepA(3), ...
                           kepA(4), kepA(5), kepA(6), mu_sun);
    
    TOF1_est = (t_GA - t_dep1_est) * 86400;
    TOF2_est = (t_arr2_est - t_GA) * 86400;
    
    [~,~,~,~,vi1,vf1,~,~] = lambertMR(rM, rE, TOF1_est, mu_sun);
    [~,~,~,~,vi2,vf2,~,~] = lambertMR(rE, rA, TOF2_est, mu_sun);
        
    dv1_est = norm(vi1' - vM);
    dv2_est = norm(vf2' - vA);
    
    % Discard worst transfers (conservative limit of 30 km/s)
    if dv1_est > 30 || dv2_est > 30
        mask(ii) = false;
        continue;
    end
        
    % Discard solutions with an excessive turning angle
    v_inf_minus = vf1' - vE;
    v_inf_plus = vi2' - vE;
        
    delta = acos(dot(v_inf_minus, v_inf_plus) / ...
                   (norm(v_inf_minus) * norm(v_inf_plus)));

    if delta > deg2rad(130)
        mask(ii) = false;
        continue;
    end
end

% Apply filters
t_GA_common = t_GA_common(mask);
K1_indices = K1_indices(mask);
K2_indices = K2_indices(mask);

fprintf('Flyby dates after filtering: %d (removed %d)\n\n', ...
    length(t_GA_common), sum(~mask)); % sum(~mask) conta quanti flyby sono stati eliminati

%% Grid search
best.DeltaV = 1e4;
best.DeltaV1 = 1e4;
best.DeltaVfb = 1e4;
best.DeltaV2 = 1e4;
best.t_dep = NaN;
best.t_GA  = NaN;
best.t_arr = NaN;
best.TOF1  = NaN;
best.TOF2  = NaN;
best.r_p = NaN;

total_iterations = length(t_GA_common);
fprintf('---- GRID SEARCH ---- \n');
fprintf('\n')

for ii = 1:length(t_GA_common)

    if mod(ii, max(1, floor(total_iterations/100))) == 0 || ii == total_iterations
        progress = (ii / total_iterations) * 100;
        fprintf('\b\b\b\b\b\b\b%6.2f%%', progress);
    end

    t_GA = t_GA_common(ii);

    t_dep1_vec = t_GA - TOF1_max : 5 : t_GA - TOF1_min;
    t_arr2_vec = t_GA + TOF2_min : 30 : t_GA + TOF2_max;

    for t_dep1 = t_dep1_vec

        tof1 = (t_GA - t_dep1) * 86400;
        if tof1 <= 0
            continue
        end

        [kepM, ~] = uplanet(t_dep1, 1);
        [rM, vM] = par2car(kepM(1), kepM(2), kepM(3), ...
                           kepM(4), kepM(5), kepM(6), mu_sun);

        [kepE, ~] = uplanet(t_GA, 3);
        [rE, vE] = par2car(kepE(1), kepE(2), kepE(3), ...
                           kepE(4), kepE(5), kepE(6), mu_sun);

        [~,~,~,~,vi1,vf1,~,~] = lambertMR(rM, rE, tof1, mu_sun);
        v_inf_minus = vf1' - vE;
        
        if norm(vi1' - vM) > 20
            continue
        end

        for t_arr2 = t_arr2_vec

            tof2 = (t_arr2 - t_GA) * 86400;
            if tof2 <= 0
                continue
            end

            [kepA,~,~] = ephAsteroids(t_arr2, 316801);
            [rA, vA] = par2car(kepA(1), kepA(2), kepA(3), ...
                               kepA(4), kepA(5), kepA(6), mu_sun);

            [~,~,~,~,vi2,vf2,~,~] = lambertMR(rE, rA, tof2, mu_sun);
            v_inf_plus = vi2' - vE;

            if norm(vf2' - vA) > 20
                continue
            end

            if norm(vi1' - vM) + norm(vf2' - vA) >= best.DeltaV
                    continue
            end

            % Turning angle
            delta = acos(dot(v_inf_minus, v_inf_plus) / ...
                       (norm(v_inf_minus) * norm(v_inf_plus)));

            if delta > deg2rad(90)
                continue
            end

            % Radius of pericenter
            fun = @(r_p) delta_powered_GA( ...
                r_p, v_inf_minus, v_inf_plus, delta, mu_earth);

            r_p_guess = 1.2 * r_p_min;
            options = optimset('TolX',1e-6,'Display','off');
            [r_p,~,exitflag] = fsolve(fun, r_p_guess, options);

            if exitflag <= 0 || r_p < r_p_min
                continue
            end

            % Semi-assi maggiori (negativi per iperboli)
            a_minus = -mu_earth / norm(v_inf_minus)^2;
            a_plus  = -mu_earth / norm(v_inf_plus)^2;

            v_p_minus = sqrt(2*mu_earth/r_p - mu_earth/a_minus);
            v_p_plus  = sqrt(2*mu_earth/r_p - mu_earth/a_plus);

            deltaV_GA = abs(v_p_plus - v_p_minus);

            deltaV = norm(vi1' - vM) + norm(vf2' - vA) + deltaV_GA;

            if deltaV < best.DeltaV
                best.DeltaV = deltaV;
                best.DeltaV1 = norm(vi1' - vM);
                best.DeltaVfb = deltaV_GA;
                best.DeltaV2 = norm(vf2' - vA);
                best.t_dep = t_dep1;
                best.t_GA  = t_GA;
                best.t_arr = t_arr2;
                best.TOF1  = tof1/86400;
                best.TOF2  = tof2/86400;
                best.r_p = r_p;
            end
        end
    end
end

fprintf('\n\nGrid search completed!\n\n');
fprintf('---- OPTIMAL SOLUTION ---- \n');
fprintf('Optimal deltaV = %.4f [km/s] \n', best.DeltaV);
fprintf('DeltaV1 (departure): %.4f [km/s] \n', best.DeltaV1);
fprintf('DeltaV_GA (flyby): %.4f [km/s] \n', best.DeltaVfb)
fprintf('DeltaV2 (arrival): %.4f [km/s] \n', best.DeltaV2)
fprintf('Departure date: %s \n', datestr(mjd20002date(best.t_dep), 'dd/mm/yyyy'));
fprintf('Flyby date: %s \n', datestr(mjd20002date(best.t_GA), 'dd/mm/yyyy'));
fprintf('Arrival date: %s \n', datestr(mjd20002date(best.t_arr), 'dd/mm/yyyy'));
fprintf('TOF1 = %.4f [days] \n', best.TOF1);
fprintf('TOF2 = %.4f [days] \n', best.TOF2);
fprintf('Total TOF = %.4f [days] \n', best.TOF1 + best.TOF2);
fprintf('Flyby radius = %.4f [km] \n\n', best.r_p);

% Find out which time window the optimal solution belongs to
best_K1 = floor((best.t_dep - t_dep1_0 + 20) / T_syn1);
best_K2 = floor((best.t_GA - t_dep2_0 + 40) / T_syn2);
fprintf('---- OPTIMAL SOLUTION WINDOW ---- \n');
fprintf('Belongs to window: K1 = %d, K2 = %d \n\n', best_K1, best_K2);

%% Refinement of the solution (fmincon)
% Initial values
x0 = [best.t_dep; best.t_GA; best.t_arr; best.r_p];

% Lower and upper bounds
lb = [best.t_dep - 10; best.t_GA - 10; best.t_arr - 10; r_p_min];
ub = [best.t_dep + 10; best.t_GA + 10; best.t_arr + 10; 5 * r_p_min];

% Options
options = optimoptions('fmincon', ...
'Display','off', ...
'Algorithm','interior-point', ...
'MaxIterations',1000, ...
'OptimalityTolerance',1e-8, ...
'StepTolerance',1e-10);

% Solve using fmincon
[x_opt, fval] = fmincon(@(x) dV_powered_GA(x, mu_sun, mu_earth), ...
                         x0, [], [], [], [], lb, ub, @time_constraints, options);
% Optimal results
opt.dep = x_opt(1);
opt.fb  = x_opt(2);
opt.arr = x_opt(3);
opt.r_p  = x_opt(4);
opt.deltaV = fval;

fprintf('\n---- REFINED OPTIMAL SOLUTION ---- \n');
fprintf('DV_min: %.4f [km/s] \n', opt.deltaV);
fprintf('Optimal Departure Date: %s \n', datestr(mjd20002date(opt.dep), 'dd/mm/yyyy'));
fprintf('Optimal Flyby Date: %s \n', datestr(mjd20002date(opt.fb), 'dd/mm/yyyy'));
fprintf('Optimal Arrival Date: %s \n', datestr(mjd20002date(opt.arr), 'dd/mm/yyyy'));
fprintf('TOF1 = %.2f [days] \n', opt.fb - opt.dep);
fprintf('TOF2 = %.2f [days] \n', opt.arr - opt.fb);
fprintf('Total TOF = %.2f [days] \n', opt.arr - opt.dep);
fprintf('Optimal r_p: %.4f [km] \n', opt.r_p);

%% Optimal parameters
opt.TOF1 = (opt.fb - opt.dep) * 86400;
opt.TOF2 = (opt.arr - opt.fb) * 86400;

% Mercury ephemeris
[kepM_opt, ~] = uplanet(opt.dep,1);
[rM_opt, vM_opt] = par2car(kepM_opt(1),kepM_opt(2),kepM_opt(3), ...
kepM_opt(4),kepM_opt(5),kepM_opt(6),mu_sun);

% Earth ephemeris
[kepE_opt, ~] = uplanet(opt.fb,3);
[rE_opt, vE_opt] = par2car(kepE_opt(1),kepE_opt(2),kepE_opt(3), ...
kepE_opt(4),kepE_opt(5),kepE_opt(6),mu_sun);

% Asteroid 316801 ephemeris
[kepA_opt, ~] = ephAsteroids(opt.arr,316801);
[rA_opt, vA_opt] = par2car(kepA_opt(1),kepA_opt(2),kepA_opt(3), ...
kepA_opt(4),kepA_opt(5),kepA_opt(6),mu_sun);

% Lambert solutions
[~,~,~,~,vi1_opt, vf1_opt,~,~] = lambertMR(rM_opt,rE_opt,opt.TOF1,mu_sun);
[~,~,~,~,vi2_opt, vf2_opt,~,~] = lambertMR(rE_opt,rA_opt,opt.TOF2,mu_sun);

% Entry and exit velocities at hyperbola
opt.v_inf_minus = vf1_opt' - vE_opt;
opt.v_inf_plus = vi2_opt' - vE_opt;

% Total turning angle
opt.delta = acos(dot(opt.v_inf_minus,opt.v_inf_plus) / ...
(norm(opt.v_inf_minus)*norm(opt.v_inf_plus)));

fprintf('\n---- FLYBY GEOMETRY ---- \n');
fprintf('Turning angle delta: %.4f [deg] \n', rad2deg(opt.delta));
fprintf('v_inf_minus: %.4f [km/s] \n', norm(opt.v_inf_minus));
fprintf('v_inf_plus: %.4f [km/s] \n', norm(opt.v_inf_plus));

% Eccentricities
opt.e_minus = 1 + (opt.r_p * norm(opt.v_inf_minus)^2) / mu_earth;
opt.e_plus = 1 + (opt.r_p * norm(opt.v_inf_plus)^2) / mu_earth;

fprintf('Eccentricity e_minus: %.4f \n', opt.e_minus);
fprintf('Eccentricity e_plus: %.4f \n', opt.e_plus);

% Turning angles
opt.delta_minus = 2 * asin(1 / opt.e_minus);
opt.delta_plus = 2 * asin(1 / opt.e_plus);
fprintf('delta_minus: %.4f [deg] \n', rad2deg(opt.delta_minus));
fprintf('delta_plus: %.4f [deg] \n', rad2deg(opt.delta_plus));

% DeltaV
opt.deltaV1 = norm(vi1_opt' - vM_opt);
opt.deltaV2 = norm(vf2_opt' - vA_opt);
opt.a_minus = -mu_earth / norm(opt.v_inf_minus)^2;
opt.a_plus = -mu_earth / norm(opt.v_inf_plus)^2;
opt.v_p_minus = sqrt(2*mu_earth/opt.r_p - mu_earth/opt.a_minus);
opt.v_p_plus = sqrt(2*mu_earth/opt.r_p - mu_earth/opt.a_plus);
opt.deltaV_GA = abs(opt.v_p_plus - opt.v_p_minus);

fprintf('\n---- DELTA-V BREAKDOWN ---- \n');
fprintf('DeltaV1 (departure): %.4f [km/s] \n', opt.deltaV1);
fprintf('DeltaV_GA (flyby): %.4f [km/s] \n', opt.deltaV_GA);
fprintf('DeltaV2 (arrival): %.4f [km/s] \n', opt.deltaV2);
fprintf('Total Delta-V: %.4f [km/s] \n\n', opt.deltaV1 + opt.deltaV_GA + opt.deltaV2);

%% Plot the transfer trajectory for this mission
% --- set up times and options ---
mjd2k1 = opt.dep; mjd2k2 = opt.fb; mjd2k3 = opt.arr;

t1 = mjd2k1 * 24 * 3600;
t2 = mjd2k2 * 24 * 3600;
t3 = mjd2k3 * 24 * 3600;
t4 = t3 + (120*24*3600);
t_full = linspace(t1, t4, 200);
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% --- get planet and asteroid heliocentric positions and velocities ---
[RM1, VM1] = get_planet_state(mjd2k1, 1, mu_sun); % mercury
[RE1, VE1] = get_planet_state(mjd2k1, 3, mu_sun); % earth
[RA1, VA1] = get_asteroid_state(mjd2k1, 316801, mu_sun); % asteroid

% --- get satellite heliocentric positions and velocities ---
RSM = RM1;
VSM = vi1_opt; % at mercury

[RSE, ~] = get_planet_state(mjd2k2, 3, mu_sun);
VSE = vi2_opt; % at earth (after GA)

[RSA, ~] = get_asteroid_state(mjd2k3, 316801, mu_sun);
VSA = vA_opt'; % at asteroid

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
plot3(R_dep(:, 1), R_dep(:, 2), R_dep(:, 3), "y--"); % during transfer

% b. earth-asteroid leg
plot3(R_GA(:, 1), R_GA(:, 2), R_GA(:, 3), "y"); % during transfer

% c. further asteroid propagation leg
plot3(R_rv(:, 1), R_rv(:, 2), R_rv(:, 3), "y"); % during transfer

% d. planet earth orbit
plot3(R_M(:, 1), R_M(:, 2), R_M(:, 3), "r"); % during transfer

% e. planet earth orbit
plot3(R_E(:, 1), R_E(:, 2), R_E(:, 3), "g"); % during transfer

% f. asteroid orbit
plot3(R_A(:, 1), R_A(:, 2), R_A(:, 3), "b"); % during transfer

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
xM = R_M(:, 1); xE = R_E(:, 1); xA = R_A(:, 1); xS = [R_dep(:, 1); R_GA(:, 1); R_rv(:, 1);];
yM = R_M(:, 2); yE = R_E(:, 2); yA = R_A(:, 2); yS = [R_dep(:, 2); R_GA(:, 2); R_rv(:, 2);];
zM = R_M(:, 3); zE = R_E(:, 3); zA = R_A(:, 3); zS = [R_dep(:, 3); R_GA(:, 3); R_rv(:, 3);];

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

%% plot flyby hyperbola in geocentric perifocal frame
% plot the two planetocentric hyperbolic arcs
figure("Name", "Planetocentric Powered Earth Fly-By"); hold on;

% --- set up ode solver conditions ---
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);
steps = 100;
tspan = linspace(0, 4800, steps);

% --- propagate and plot incoming planetocentric arc ---
y = [[opt.r_p; 0; 0] [0; opt.v_p_minus; 0;]];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_earth), -tspan, y, options);
r_minus = Y(:, 1:3) ./ R_earth;
plot3(r_minus(:, 1), r_minus(:, 2), r_minus(:, 3), "r")

% --- propagate and plot outgoing planetocentric arc ---
y = [[opt.r_p; 0; 0] [0; opt.v_p_plus; 0;]];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_earth), tspan, y, options);
r_plus = Y(:, 1:3) ./ R_earth;
plot3(r_plus(:, 1), r_plus(:, 2), r_plus(:, 3), "g")

% --- plot earth surface texture ---
earth_img = imread("EarthTexture.jpg");
[x_earth, y_earth, z_earth] = sphere(50);
surface(x_earth, y_earth, -z_earth, "FaceColor", ...
    "texturemap", "CData", earth_img, "EdgeColor", "none")

% --- finish up plot properties ---
legend("Incoming trajectory", "Outgoing trajectory", "", "Sun position")
xlabel("X [r_E_a_r_t_h]"); ylabel("Y [r_E_a_r_t_h]"); zlabel("Z [r_E_a_r_t_h]");
title("Powered Earth Fly-by Trajectory");
grid on; axis equal; hold off;

%% === part 5 (additional) ===
% i want to plot the heliocentric trajectory too
figure("Name", "Heliocentric Powered Earth Fly-By"); hold on;
steps = 1000;
tspan_helio = linspace(0, 60*60*24*30, steps);

% --- plot position of earth during fly-by ---
y = [RE1 VE1];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), -tspan_helio, y, options);
R_E_minus = Y(:, 1:3);
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), 9*tspan_helio, y, options);
R_E_plus = Y(:, 1:3);
plot3(R_E_minus(:, 1), R_E_minus(:, 2), R_E_minus(:, 3), "g");
plot3(R_E_plus(:, 1), R_E_plus(:, 2), R_E_plus(:, 3), "g");
scatter3(RSE(1), RSE(2), RSE(3), "filled", "g");

% --- plot earth surface texture ---
earth_img = imread("EarthTexture.jpg");
[x_earth, y_earth, z_earth] = sphere(50);
surface(x_earth, y_earth, -z_earth, "FaceColor", ...
    "texturemap", "CData", earth_img, "EdgeColor", "none")

% --- propagate and plot incoming heliocentric arc ---
% y = [RSE vf1_opt];
% [~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), -5*tspan_helio, y, options);
% R_minus = Y(:, 1:3);
% plot3(R_minus(:, 1), R_minus(:, 2), R_minus(:, 3), "r")
R_dep = R_dep .* AU;
plot3(R_dep(:, 1), R_dep(:, 2), R_dep(:, 3), "y--"); % incoming trajectory

% --- propagate and plot outgoing heliocentric arc ---
y = [RSE VSE];
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_sun), 5*tspan_helio, y, options);
R_plus = Y(:, 1:3);
plot3(R_plus(:, 1), R_plus(:, 2), R_plus(:, 3), "y") % outgoing trajectory

% --- plot sun ---
scatter3(0, 0, 0, "filled", "y");

% --- finish up plot properties ---
legend("", "", "Earth", "Incoming Trajectory", "Outgoing Trajectory", "Sun")
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Powered Earth Fly-by Trajectory");
grid on; axis equal; hold off;

















%% SELF NOTES SELF NOTES SELF NOTES SELF NOTES SELF NOTES SELF NOTES
% porkchop plot ME and EA refined solution (using fmincon)
% FATTO plot 3D total solution 
% FATTO FILMATO CLAMOROSO (da lasciare solo su matlab, ma da mettere
% assolutamente)
% plot hyperbola
% commentare e descrivere come vogliono loro tutte le funzioni:
% - dV_powered_GA
% - delta_powered_GA
% - time_constraints
% modificare par2car ed eventualmente altro (Terra3D, ecc)
% report (max 6 pages):
% - introduction
% - circular and coplanar approx --> T_syn and TOF limits
% - porkchop plots first T_syn
% - design strategy (to find t_GA_common)
% - grid search
% - refinement and final results
% - plots and final remarks
