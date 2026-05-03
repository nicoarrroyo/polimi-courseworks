
clc; clear; close all;
fprintf('=================================================================\n');
fprintf('  ADCS REVERSE SIZING  –  MARS EXPRESS (MEX)\n');
fprintf('=================================================================\n\n');

%% ========================================================================
%  1. MISSION & ENVIRONMENT PARAMETERS
% =========================================================================

% --- Spacecraft parameters -----------------------------------------------
dim       = [1.5, 1.8, 1.0];           % Dimensions [m]  (x, y, z)
Ix        = 373.08;                    % Moment of inertia – X axis [kg m²]
Iy        = 1012.37;                   % Moment of inertia – Y axis [kg m²]
Iz        = 1088.31;                   % Moment of inertia – Z axis [kg m²]
I_all     = [Ix, Iy, Iz];
I_min     = min(I_all);                % Minimum inertia [kg m²]
I_max     = max(I_all);                % Maximum inertia [kg m²]

A_cross   = 11.9;                      % Cross-sectional area [m²]
d_com     = 0.1;                       % CoM displacement (csp – cg) [m]
q_refl    = 0.3;                       % Reflectivity coefficient [-]

T_rw_max  = 0.075;                     % RW max torque [Nm]
h_rw_max  = 12.0;                      % RW max momentum [Nms]

F_thr     = 10;                        % Thruster force [N]
L_thr     = 0.69;                      % Thruster moment arm [m]
n_thr     = 2;                         % Number of thrusters in couple [-]
Isp       = 292;                       % Specific impulse [s]
g0        = 9.81;                      % Standard gravity [m/s²]

rate_max  = 0.5 * (pi/180);            % Max angular rate [rad/s]
theta_max = 180 * (pi/180);            % Worst-case slew angle [rad]

% --- Orbital / environment parameters ------------------------------------
a         = 9354.1e3;                  % Semi-major axis [m]
T_orb     = 2.75e4;                    % Orbit period [s]
mu        = 4.28e13;                   % Mars gravitational constant [m³/s²]
R_Mars    = 3396.2e3;                  % Mars equatorial radius [m]
F_sun     = 590;                       % Solar constant at Mars [W/m²]
c_light   = 3e8;                       % Speed of light [m/s]

% Orbital radius (approximated as semi-major axis – conservative)
R_orb     = a;

fprintf('--- Spacecraft Parameters ---\n');
fprintf('  Inertia:   Ix = %.2f | Iy = %.2f | Iz = %.2f  kg m²\n', Ix, Iy, Iz);
fprintf('  I_min = %.2f kg m²   I_max = %.2f kg m²\n', I_min, I_max);
fprintf('  RW max torque    = %.3f Nm\n', T_rw_max);
fprintf('  RW max momentum  = %.1f Nms\n', h_rw_max);
fprintf('  Thruster force   = %.1f N  |  arm = %.2f m  |  n = %d\n\n', F_thr, L_thr, n_thr);

%% ========================================================================
%  2. DISTURBANCE ENVIRONMENT EVALUATION
% =========================================================================

fprintf('=================================================================\n');
fprintf('  2. DISTURBANCE ENVIRONMENT\n');
fprintf('=================================================================\n\n');

% -------------------------------------------------------------------------
%  2.1  Gravity Gradient Torque
% -------------------------------------------------------------------------
% Worst-case: payload (aligned with I_max axis) pointed 22° off-nadir
theta_gg = 22 * (pi/180);             % Off-nadir angle [rad]

T_gg = (3*mu) / (2*R_orb^3) * (I_max - I_min) * sin(2*theta_gg);

fprintf('--- 2.1  Gravity Gradient Torque ---\n');
fprintf('  Off-nadir angle (worst case) : theta = %.1f°\n', theta_gg*(180/pi));
fprintf('  T_gg = 3µ/(2R³) · (I_max – I_min) · sin(2θ)\n');
fprintf('  T_gg = %.4e  Nm\n\n', T_gg);

% -------------------------------------------------------------------------
%  2.2  Solar Radiation Pressure Torque
% -------------------------------------------------------------------------
% Worst-case: Sun incidence angle i = 0°  →  cos(0) = 1
i_sun = 0;                             % Sun incidence angle [rad]

T_srp = (F_sun / c_light) * A_cross * (1 + q_refl) * cos(i_sun) * d_com;

fprintf('--- 2.2  Solar Radiation Pressure Torque ---\n');
fprintf('  Sun incidence angle (worst case) : i = %.1f°\n', i_sun*(180/pi));
fprintf('  T_SRP = (Fs/c) · A · (1+q) · cos(i) · (csp–cg)\n');
fprintf('  T_SRP = %.4e  Nm\n\n', T_srp);

fprintf('  Ratio T_gg / T_SRP = %.1f  (SRP ~one order of magnitude smaller)\n\n', T_gg/T_srp);

%% ========================================================================
%  3. REACTION WHEEL SIZING
% =========================================================================

fprintf('=================================================================\n');
fprintf('  3. REACTION WHEEL SIZING\n');
fprintf('=================================================================\n\n');

% -------------------------------------------------------------------------
%  3.1  Momentum accumulation due to disturbances
% -------------------------------------------------------------------------
% Only CONSTANT disturbances drive secular buildup (cyclic disturbances
% average to zero over one orbit).
% For a nadir-pointing spacecraft: T_gg is constant, T_SRP is cyclic.
% Conservative factor of 2 is applied to T_gg.

T_dis = 2 * T_gg;                     % Total constant disturbance torque [Nm]

% Momentum accumulated per orbit
h_orbit = T_dis * T_orb;              % [Nms]

% Number of orbits before RW saturation
N_orb = h_rw_max / h_orbit;

fprintf('--- 3.1  Disturbance Momentum Accumulation ---\n');
fprintf('  T_dis = 2 × T_gg = %.4e  Nm\n', T_dis);
fprintf('  h_orbit = T_dis × T_orb = %.4f  Nms\n', h_orbit);
fprintf('  N_orb (before saturation) = h_max / h_orbit = %.2f  orbits\n\n', N_orb);

% Check cyclic SRP peak (should be negligible)
omega_orb    = 2*pi / T_orb;
h_cyclic_max = (T_srp / omega_orb) * 2;   % integral of sin over half-period × 2
fprintf('  Cyclic SRP peak momentum = %.4f  Nms  (≪ h_max = %.1f Nms ✓)\n\n', h_cyclic_max, h_rw_max);

% -------------------------------------------------------------------------
%  Plot: Momentum accumulation over 5.6 orbits
% -------------------------------------------------------------------------
t_plot   = linspace(0, N_orb * T_orb, 5000);
h_sec    = T_dis * t_plot; % Secular (GG)
h_cyc    = (T_srp / omega_orb) * (1 - cos(omega_orb*t_plot)); % Cyclic (SRP)
h_total  = h_sec + h_cyc;

figure('Name','Momentum Accumulation','Color','w','Position',[100 100 800 420]);
hold on; grid on;
plot(t_plot/T_orb, h_sec,   '--b', 'LineWidth', 1.5, 'DisplayName','Secular Buildup (GG)');
plot(t_plot/T_orb, h_cyc,   '--r', 'LineWidth', 1.2, 'DisplayName','Cyclic Variation (SRP)');
plot(t_plot/T_orb, h_total, '-k',  'LineWidth', 2.0, 'DisplayName','Total Stored Momentum');
yline(h_rw_max, '--m', 'LineWidth', 1.5, 'DisplayName', 'h_{max}');
xlabel('Time [Orbits]', 'FontSize', 12);
ylabel('Momentum [Nms]', 'FontSize', 12);
title('Reaction Wheel Momentum Loading', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location','northwest', 'FontSize', 10);
xlim([0, N_orb]); ylim([0, h_rw_max*1.1]);

% -------------------------------------------------------------------------
%  3.2  Slew Maneuver Sizing
% -------------------------------------------------------------------------
fprintf('--- 3.2  Slew Maneuver ---\n');
fprintf('  Worst-case slew angle : θ_slew = %.0f°\n', theta_max*(180/pi));

% Minimum slew time using maximum torque
t_slew_min = sqrt(4 * theta_max * I_max / T_rw_max);
fprintf('\n  [Case A] Maximum torque applied:\n');
fprintf('  t_slew_min = sqrt(4·θ·I_max / T_max) = %.2f  s\n', t_slew_min);

% Peak angular rate and momentum at max torque
rate_slew_max_A = (T_rw_max / I_max) * (t_slew_min / 2);
h_slew_max_A    = I_max * rate_slew_max_A;
fprintf('  θ̇_max = T_max/I_max · (t/2) = %.4f  rad/s  (%.4f °/s)\n', ...
        rate_slew_max_A, rate_slew_max_A*(180/pi));
fprintf('  h_slew_peak = I_max · θ̇_max = %.2f  Nms\n', h_slew_max_A);
fprintf('  ⚠  Angular rate limit EXCEEDED (%.2f°/s > %.1f°/s) – torque must be reduced.\n', ...
        rate_slew_max_A*(180/pi), rate_max*(180/pi));

% Reduced torque to comply with angular rate constraint
% θ̇_max = T/I_max · t_slew/2  and  θ_slew = T/I_max · (t_slew/2)²
% → T = I_max · θ̇_max² / θ_slew (from inversion of the parabolic profile)
T_reduced  = I_max * rate_max^2 / theta_max;
t_slew_red = 2 * rate_max * I_max / T_reduced;
h_slew_max_B = I_max * rate_max;

fprintf('\n  [Case B] Reduced torque to satisfy θ̇ ≤ %.1f°/s:\n', rate_max*(180/pi));
fprintf('  T_reduced  = I_max · θ̇_max² / θ_slew = %.4f  Nm\n', T_reduced);
fprintf('  t_slew     = 2 · θ̇_max · I_max / T   = %.1f  s\n', t_slew_red);
fprintf('  h_slew_peak = I_max · θ̇_max           = %.2f  Nms\n\n', h_slew_max_B);

% -------------------------------------------------------------------------
%  Plot: Angular rate profiles for both cases
% -------------------------------------------------------------------------
figure('Name','Slew Angular Rate Profiles','Color','w','Position',[100 550 800 420]);
hold on; grid on;

% Case A – max torque (trapezoidal, symmetric)
t_A = [0, t_slew_min/2, t_slew_min];
r_A = [0, rate_slew_max_A, 0] * (180/pi);
plot(t_A, r_A, '-r', 'LineWidth', 2, 'DisplayName', sprintf('Max torque (T=%.3f Nm)', T_rw_max));

% Case B – reduced torque
t_B = [0, t_slew_red/2, t_slew_red];
r_B = [0, rate_max, 0] * (180/pi);
plot(t_B, r_B, '-b', 'LineWidth', 2, 'DisplayName', sprintf('Reduced torque (T=%.3f Nm)', T_reduced));

yline(rate_max*(180/pi), '--k', 'LineWidth', 1.2, 'DisplayName', 'Rate limit 0.5°/s');
xlabel('Time [s]', 'FontSize', 12);
ylabel('Angular Rate [°/s]', 'FontSize', 12);
title('Slew Maneuver – Angular Rate Profiles (Reaction Wheels)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'north', 'FontSize', 10);

%% ========================================================================
%  4. THRUSTER SIZING – SLEW MANEUVER
% =========================================================================

fprintf('=================================================================\n');
fprintf('  4. THRUSTER SIZING – SLEW MANEUVER\n');
fprintf('=================================================================\n\n');

% Maneuver: acceleration + coasting + braking
% Constraint: peak angular rate = rate_max

% Acceleration phase duration
t_acc = (rate_max * I_max) / (n_thr * F_thr * L_thr);

% Angle spanned during acceleration (= same for braking)
theta_acc = 0.5 * (n_thr * F_thr * L_thr / I_max) * t_acc^2;

% Coasting angle and time
theta_coast = theta_max - 2*theta_acc;
t_coast     = theta_coast / rate_max;

% Total slew time
t_slew_thr  = t_acc + t_coast + t_acc;

fprintf('  Acceleration/braking time  t_acc   = %.4f  s\n', t_acc);
fprintf('  Angle per phase            θ_acc   = %.6f  rad  (%.4f°)\n', theta_acc, theta_acc*(180/pi));
fprintf('  Coasting angle             θ_coast = %.4f  rad  (%.4f°)\n', theta_coast, theta_coast*(180/pi));
fprintf('  Coasting time              t_coast = %.2f  s\n', t_coast);
fprintf('  Total slew time            t_slew  = %.2f  s\n\n', t_slew_thr);

% -------------------------------------------------------------------------
%  Plot: Thruster slew angular rate profile
% -------------------------------------------------------------------------
t_thr_full = [0, t_acc, t_acc+t_coast, t_slew_thr];
r_thr      = [0, rate_max, rate_max, 0] * (180/pi);

figure('Name','Thruster Slew Profile','Color','w','Position',[100 100 800 400]);
hold on; grid on;
plot(t_thr_full, r_thr, '-g', 'LineWidth', 2.5);
yline(rate_max*(180/pi), '--k', 'LineWidth', 1.2, 'DisplayName','Rate limit 0.5°/s');
xlabel('Time [s]', 'FontSize', 12);
ylabel('Angular Rate [°/s]', 'FontSize', 12);
title('Slew Maneuver – Angular Rate Profile (Thrusters)', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Thruster slew profile','Rate limit 0.5°/s'}, 'Location','northeast','FontSize',10);

%% ========================================================================
%  5. PROPELLANT MASS SIZING – WHEEL DESATURATION
% =========================================================================

fprintf('=================================================================\n');
fprintf('  5. PROPELLANT MASS SIZING – WHEEL DESATURATION\n');
fprintf('=================================================================\n\n');

% Theoretical minimum desaturation time (momentum transfer at full thrust)
t_desat_min = h_rw_max / (n_thr * L_thr * F_thr);
fprintf('  Theoretical minimum desaturation time = %.2f  s\n', t_desat_min);
fprintf('  (Requires RW torque >> T_max → physically unfeasible)\n\n');

% Realistic desaturation: driven by RW max torque (use T_max/2 conservatively)
% Time to dump h_max at T_max/2 from each side (momentum transferred by thrusters)
t_desat = 2 * h_rw_max / T_rw_max;

% Required force per thruster (couple: n = 2)
F_desat = h_rw_max / (n_thr * L_thr * t_desat);

fprintf('  Realistic desaturation (limited by T_max/2):\n');
fprintf('  t_desat  = 2·h_max / T_max = %.1f  s\n', t_desat);
fprintf('  F_desat  = h_max / (n·L·t) = %.4f  N\n\n', F_desat);

% Number of desaturation manoeuvres over 2-year mission
T_mission = 2 * 365.25 * 24 * 3600;   % [s]
N_desat   = (T_mission / T_orb) / N_orb;
fprintf('  Mission duration    = 2 years = %.2e  s\n', T_mission);
fprintf('  Total orbits        = %.0f\n', T_mission / T_orb);
fprintf('  Desaturation events = %.0f\n\n', N_desat);

% Propellant mass (Tsiolkovsky / rocket equation form for impulsive burns)
% M_prop = I_total / (Isp · g0),   I_total = N_desat · n · t_desat · F_desat
I_total  = N_desat * n_thr * t_desat * F_desat;
M_prop   = I_total / (Isp * g0);

fprintf('  Total impulse for desaturation:\n');
fprintf('  I_total = N_desat × n × t_desat × F_desat = %.2f  Ns\n', I_total);
fprintf('  Propellant mass M_prop = I_total / (Isp·g0) = %.2f  kg\n\n', M_prop);
fprintf('  ⚠  Propulsion subsystem margins shall be applied on top of this value.\n\n');

%% ========================================================================
%  6. SUMMARY TABLE
% =========================================================================

fprintf('=================================================================\n');
fprintf('  6. RESULTS SUMMARY\n');
fprintf('=================================================================\n\n');

fprintf('  %-45s  %12.4e  Nm\n',  'Gravity Gradient Torque T_gg',          T_gg);
fprintf('  %-45s  %12.4e  Nm\n',  'SRP Torque T_srp',                       T_srp);
fprintf('  %-45s  %12.4e  Nm\n',  'Total constant disturbance T_dis',        T_dis);
fprintf('  %-45s  %12.4f   Nms\n','Momentum per orbit h_orbit',              h_orbit);
fprintf('  %-45s  %12.2f   orb\n','Orbits before RW saturation N_orb',       N_orb);
fprintf('\n');
fprintf('  %-45s  %12.2f   s\n',  'Min slew time (max torque) t_slew_min',   t_slew_min);
fprintf('  %-45s  %12.4f   °/s\n','Peak rate at max torque',                 rate_slew_max_A*(180/pi));
fprintf('  %-45s  %12.4f   Nm\n', 'Reduced torque (rate-constrained)',        T_reduced);
fprintf('  %-45s  %12.1f   s\n',  'Slew time (rate-constrained)',             t_slew_red);
fprintf('  %-45s  %12.2f   Nms\n','Peak RW momentum during slew',             h_slew_max_B);
fprintf('\n');
fprintf('  %-45s  %12.4f   s\n',  'Thruster slew t_acc',                     t_acc);
fprintf('  %-45s  %12.2f   s\n',  'Thruster slew total time',                 t_slew_thr);
fprintf('\n');
fprintf('  %-45s  %12.1f   s\n',  'Realistic desaturation duration',          t_desat);
fprintf('  %-45s  %12.4f   N\n',  'Thruster force for desaturation',          F_desat);
fprintf('  %-45s  %12.0f\n',      'Total desaturation events (2-yr mission)', N_desat);
fprintf('  %-45s  %12.2f   kg\n', 'Required propellant mass',                 M_prop);
fprintf('\n=================================================================\n');
fprintf('  END OF ADCS SIZING\n');
fprintf('=================================================================\n');
