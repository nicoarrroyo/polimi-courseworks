
clear; close all; clc;
fprintf('===========================================================\n');
fprintf('  ADCS REVERSE SIZING  –  MARS EXPRESS (MEX)\n');
fprintf('===========================================================\n\n');

%% 1. KNOWN PARAMETERS (Masses in kg, Dimensions in meters)
% MIRAS Arms
m_arm           = 90.5;
m_arms_total    = m_arm * 3;
arm_L           = 4.6;

% MIRAS Hub
m_hub           = 83.5;
hub_d           = 1.3;
hub_r           = hub_d / 2;    % 0.65 m
hub_lz          = 1.2;

% Solar Panels
m_solar         = 49.92;
solar_lx        = 0.01;         % Approximated as near-zero thickness
solar_ly        = 4.0;
solar_lz        = 1.5;

% Proteus Bus
m_bus           = 225.08;
bus_lx          = 0.954;
bus_ly          = 0.954;
bus_lz          = 1.0;

%% 2. MISSING VARIABLES
% Note: All Z-coordinates should be relative to the same datum (e.g., bottom of the bus = 0)

z_arms_com  = 0; % Z-coordinate for the Arms Center of Mass (ASSUMING 0 thickness)
z_hub_com   = hub_lz / 2; % Z-coordinate Hub Center of Mass
z_solar_com = hub_lz + (solar_lz / 2); % Z-coordinate Solar Panel Center of Mass
z_bus_com   = hub_lz + solar_lz + (bus_lz / 2); % Z-coordinate Bus Center of Mass

r_arm_start = 0.65; % Radial offset where arms attach

%% 3. CALCULATE LOCAL INERTIA TENSORS (About their own COMs)

% Bus (Rectangular Prism)
I_bus_local = diag([ ...
    (m_bus/12) * (bus_ly^2 + bus_lz^2), ... % Ixx
    (m_bus/12) * (bus_lx^2 + bus_lz^2), ... % Iyy
    (m_bus/12) * (bus_lx^2 + bus_ly^2)  ... % Izz
]);

% Solar Panels (Rectangular Prism)
I_solar_local = diag([ ...
    (m_solar/12) * (solar_ly^2 + solar_lz^2), ...
    (m_solar/12) * (solar_lx^2 + solar_lz^2), ...
    (m_solar/12) * (solar_lx^2 + solar_ly^2)  ...
]);

% MIRAS Hub (Solid Cylinder Approximation)
I_hub_local = diag([ ...
    (m_hub/12) * (3*hub_r^2 + hub_lz^2), ...
    (m_hub/12) * (3*hub_r^2 + hub_lz^2), ...
    (m_hub/2)  * (hub_r^2) ...
]);

% MIRAS Arms (3 Rods offset by r_arm_start, spread 120 deg)
% Inertia of a single rod about the central Z-axis (evaluating integral of x^2 dm from r0 to r0+L)
I_single_arm = (m_arm / 3) * (arm_L^2 + 3*r_arm_start*arm_L + 3*r_arm_start^2);

% Summing 3 arms at 0, 120, and 240 degrees perfectly cancels cross-products.
I_arms_local = diag([ ...
    1.5 * I_single_arm, ... % Ixx
    1.5 * I_single_arm, ... % Iyy
    3.0 * I_single_arm  ... % Izz
]);

%% 4. CALCULATE GLOBAL CENTER OF MASS (Z-Axis)
total_mass = m_bus + m_solar + m_hub + m_arms_total;

z_cg_global = (m_bus * z_bus_com + ...
               m_solar * z_solar_com + ...
               m_hub * z_hub_com + ...
               m_arms_total * z_arms_com) / total_mass;

%% 5. APPLY PARALLEL AXIS THEOREM (Shift to Global COM)
% Function to apply Z-shift to an inertia tensor
% Shifting along Z only affects Ixx and Iyy. Izz remains unchanged.
shift_tensor = @(mass, dz) diag([mass*dz^2, mass*dz^2, 0]);

dz_bus   = z_bus_com - z_cg_global;
dz_solar = z_solar_com - z_cg_global;
dz_hub   = z_hub_com - z_cg_global;
dz_arms  = z_arms_com - z_cg_global;

I_bus_shifted   = I_bus_local   + shift_tensor(m_bus, dz_bus);
I_solar_shifted = I_solar_local + shift_tensor(m_solar, dz_solar);
I_hub_shifted   = I_hub_local   + shift_tensor(m_hub, dz_hub);
I_arms_shifted  = I_arms_local  + shift_tensor(m_arms_total, dz_arms);

I_total = I_bus_shifted + I_solar_shifted + I_hub_shifted + I_arms_shifted;

%% 6. DISPLAY RESULTS
fprintf('========================================\n');
fprintf(' SMOS SATELLITE INERTIA CALCULATION\n');
fprintf('========================================\n');
fprintf(' Total System Mass:  %.2f kg\n', total_mass);
fprintf(' Global Center of Mass (Z): %.4f m\n\n', z_cg_global);

fprintf(' Total Inertia Tensor (kg*m^2) at Global COM:\n');
fprintf(' Ixx: %10.2f\n', I_total(1,1));
fprintf(' Iyy: %10.2f\n', I_total(2,2));
fprintf(' Izz: %10.2f\n', I_total(3,3));
fprintf('\nNote: Products of inertia (Ixy, Ixz, Iyz) are assumed 0 \ndue to geometric symmetry and assumptions.\n');
fprintf('========================================\n');

%% ========================================================================
%  1. MISSION & ENVIRONMENT PARAMETERS
% =========================================================================

% --- Spacecraft parameters -----------------------------------------------
dim       = [1.5, 1.8, 1.0];    % Dimensions [m]  (x, y, z)
Ix        = I_total(1,1);       % Moment of inertia – X axis [kg m²]
Iy        = I_total(2,2);       % Moment of inertia – Y axis [kg m²]
Iz        = I_total(3,3);       % Moment of inertia – Z axis [kg m²]
I_all     = [Ix, Iy, Iz];
I_min     = min(I_all);         % Minimum inertia [kg m²]
I_max     = max(I_all);         % Maximum inertia [kg m²]

A_cross   = 11.9;               % Cross-sectional area [m²]
d_com     = 0.1;                % CoM displacement (csp – cg) [m]
q_refl    = 0.3;                % Reflectivity coefficient [-]

T_rw_max  = 0.075;              % RW max torque [Nm]
h_rw_max  = 12.0;               % RW max momentum [Nms]

F_thr     = 10;                 % Thruster force [N]
L_thr     = 0.69;               % Thruster moment arm [m]
n_thr     = 2;                  % Number of thrusters in couple [-]
Isp       = 292;                % Specific impulse [s]
g0        = 9.81;               % Standard gravity [m/s²]

rate_max  = 0.5 * (pi/180);     % Max angular rate [rad/s]
theta_max = 180 * (pi/180);     % Worst-case slew angle [rad]

% --- Orbital / environment parameters ------------------------------------
a         = 9354.1e3;           % Semi-major axis [m]
T_orb     = 2.75e4;             % Orbit period [s]
mu        = 398600.4;           % Mars gravitational constant [m³/s²]
R_Mars    = 3396.2e3;           % Mars equatorial radius [m]
F_sun     = 590;                % Solar constant at Mars [W/m²]
c_light   = 3e8;                % Speed of light [m/s]

% Orbital radius (approximated as semi-major axis – conservative)
R_orb     = a;

fprintf('--- Spacecraft Parameters ---\n');
fprintf(' Inertia:               Ix = %.2f kg m^2\n', Ix);
fprintf(' Inertia:               Iy = %.2f kg m^2\n', Iy);
fprintf(' Inertia:               Iz = %.2f kg m^2\n', Iz);
fprintf(' Inertia:            I_min = %.2f kg m^2\n', I_min);
fprintf(' Inertia:            I_max = %.2f kg m^2\n', I_max);
fprintf(' RW:            max torque = %.3f Nm\n', T_rw_max);
fprintf(' RW:          max momentum = %.1f Nms\n', h_rw_max);
fprintf(' Thruster:           force = %.1f N\n', F_thr);
fprintf(' Thruster:             arm = %.2f m\n', L_thr);
fprintf(' Thruster:               n = %d\n\n', n_thr);

%% ========================================================================
%  2. DISTURBANCE ENVIRONMENT EVALUATION
% =========================================================================

fprintf('===========================================================\n');
fprintf('  2. DISTURBANCE ENVIRONMENT\n');
fprintf('===========================================================\n\n');

% -------------------------------------------------------------------------
%  2.1  Gravity Gradient Torque
% -------------------------------------------------------------------------
% Worst-case: payload (aligned with I_max axis) pointed 22 deg off-nadir
theta_gg = 22 * (pi/180);   % Off-nadir angle [rad]

T_gg = (3*mu) / (2*R_orb^3) * (I_max - I_min) * sin(2*theta_gg);

fprintf('--- 2.1  Gravity Gradient Torque ---\n');
fprintf(' Worst-case off-nadir angle : theta = %.1f deg\n', theta_gg*(180/pi));
fprintf(' T_gg = %.4e Nm\n\n', T_gg);

% -------------------------------------------------------------------------
%  2.2  Solar Radiation Pressure Torque
% -------------------------------------------------------------------------
% Worst-case: Sun incidence angle i = 0 deg -> cos(0) = 1
i_sun = 0;                  % Sun incidence angle [rad]

T_srp = (F_sun / c_light) * A_cross * (1 + q_refl) * cos(i_sun) * d_com;

fprintf('--- 2.2  Solar Radiation Pressure Torque ---\n');
fprintf(' Worst case Sun incidence angle : i = %.1f deg\n', i_sun*(180/pi));
fprintf(' T_SRP = %.4e Nm\n\n', T_srp);

fprintf(' Ratio T_gg / T_SRP = %.1f\n\n', T_gg/T_srp);

%% ========================================================================
%  3. REACTION WHEEL SIZING
% =========================================================================

fprintf('===========================================================\n');
fprintf('  3. REACTION WHEEL SIZING\n');
fprintf('===========================================================\n\n');

% -------------------------------------------------------------------------
%  3.1  Momentum accumulation due to disturbances
% -------------------------------------------------------------------------
% Only CONSTANT disturbances drive secular buildup (cyclic disturbances
% average to zero over one orbit).
% For a nadir-pointing spacecraft: T_gg is constant, T_SRP is cyclic.
% Conservative factor of 2 is applied to T_gg.

T_dis = 2 * T_gg;               % Total constant disturbance torque [Nm]

h_orbit = T_dis * T_orb;        % Momentum accumulated per orbit [Nms]

% Number of orbits before RW saturation
N_orb = h_rw_max / h_orbit;

fprintf('--- 3.1  Disturbance Momentum Accumulation ---\n');
fprintf(' T_dis = 2 × T_gg = %.4e Nm\n', T_dis);
fprintf(' h_orbit = T_dis × T_orb = %.4f Nms\n', h_orbit);
fprintf(' N_orb (before saturation) = h_max / h_orbit = %.2f  orbits\n\n', N_orb);

% Check cyclic SRP peak (should be negligible)
omega_orb    = 2*pi / T_orb;
h_cyclic_max = (T_srp / omega_orb) * 2;   % integral of sin over half-period × 2
fprintf(' Cyclic SRP peak momentum = %.4f  Nms  (≪ h_max = %.1f Nms ✓)\n\n', h_cyclic_max, h_rw_max);

% -------------------------------------------------------------------------
%  Plot: Momentum accumulation over 5.6 orbits
% -------------------------------------------------------------------------
t_plot   = linspace(0, N_orb * T_orb, 5000);
h_sec    = T_dis * t_plot; % Secular (GG)
h_cyc    = (T_srp / omega_orb) * (1 - cos(omega_orb*t_plot)); % Cyclic (SRP)
h_total  = h_sec + h_cyc;

% figure('Name','Momentum Accumulation','Color','w','Position',[100 100 800 420]);
% hold on; grid on;
% plot(t_plot/T_orb, h_sec,   '--b', 'LineWidth', 1.5, 'DisplayName','Secular Buildup (GG)');
% plot(t_plot/T_orb, h_cyc,   '--r', 'LineWidth', 1.2, 'DisplayName','Cyclic Variation (SRP)');
% plot(t_plot/T_orb, h_total, '-k',  'LineWidth', 2.0, 'DisplayName','Total Stored Momentum');
% yline(h_rw_max, '--m', 'LineWidth', 1.5, 'DisplayName', 'h_{max}');
% xlabel('Time [Orbits]', 'FontSize', 12);
% ylabel('Momentum [Nms]', 'FontSize', 12);
% title('Reaction Wheel Momentum Loading', 'FontSize', 13, 'FontWeight', 'bold');
% legend('Location','northwest', 'FontSize', 10);
% xlim([0, N_orb]); ylim([0, h_rw_max*1.1]);

% -------------------------------------------------------------------------
%  3.2  Slew Maneuver Sizing
% -------------------------------------------------------------------------
fprintf('--- 3.2  Slew Maneuver ---\n');
fprintf(' Worst-case slew angle : θ_slew = %.0f°\n', theta_max*(180/pi));

% Minimum slew time using maximum torque
t_slew_min = sqrt(4 * theta_max * I_max / T_rw_max);
fprintf('\n  [Case A] Maximum torque applied:\n');
fprintf(' t_slew_min = sqrt(4·θ·I_max / T_max) = %.2f  s\n', t_slew_min);

% Peak angular rate and momentum at max torque
rate_slew_max_A = (T_rw_max / I_max) * (t_slew_min / 2);
h_slew_max_A    = I_max * rate_slew_max_A;
fprintf(' θ̇_max = T_max/I_max · (t/2) = %.4f  rad/s  (%.4f °/s)\n', ...
        rate_slew_max_A, rate_slew_max_A*(180/pi));
fprintf(' h_slew_peak = I_max · θ̇_max = %.2f  Nms\n', h_slew_max_A);
fprintf(' !! Angular rate limit EXCEEDED (%.2f°/s > %.1f°/s) – torque must be reduced.\n', ...
        rate_slew_max_A*(180/pi), rate_max*(180/pi));

% Reduced torque to comply with angular rate constraint
% θ̇_max = T/I_max · t_slew/2  and  θ_slew = T/I_max · (t_slew/2)²
% → T = I_max · θ̇_max² / θ_slew (from inversion of the parabolic profile)
T_reduced  = I_max * rate_max^2 / theta_max;
t_slew_red = 2 * rate_max * I_max / T_reduced;
h_slew_max_B = I_max * rate_max;

fprintf('\n  [Case B] Reduced torque to satisfy θ̇ ≤ %.1f°/s:\n', rate_max*(180/pi));
fprintf(' T_reduced  = I_max · θ̇_max² / θ_slew = %.4f  Nm\n', T_reduced);
fprintf(' t_slew     = 2 · θ̇_max · I_max / T   = %.1f  s\n', t_slew_red);
fprintf(' h_slew_peak = I_max · θ̇_max           = %.2f  Nms\n\n', h_slew_max_B);

% -------------------------------------------------------------------------
%  Plot: Angular rate profiles for both cases
% -------------------------------------------------------------------------
% figure('Name','Slew Angular Rate Profiles','Color','w','Position',[100 550 800 420]);
% hold on; grid on;

% Case A – max torque (trapezoidal, symmetric)
% t_A = [0, t_slew_min/2, t_slew_min];
% r_A = [0, rate_slew_max_A, 0] * (180/pi);
% plot(t_A, r_A, '-r', 'LineWidth', 2, 'DisplayName', sprintf('Max torque (T=%.3f Nm)', T_rw_max));

% Case B – reduced torque
% t_B = [0, t_slew_red/2, t_slew_red];
% r_B = [0, rate_max, 0] * (180/pi);
% plot(t_B, r_B, '-b', 'LineWidth', 2, 'DisplayName', sprintf('Reduced torque (T=%.3f Nm)', T_reduced));
% 
% yline(rate_max*(180/pi), '--k', 'LineWidth', 1.2, 'DisplayName', 'Rate limit 0.5°/s');
% xlabel('Time [s]', 'FontSize', 12);
% ylabel('Angular Rate [°/s]', 'FontSize', 12);
% title('Slew Maneuver – Angular Rate Profiles (Reaction Wheels)', 'FontSize', 13, 'FontWeight', 'bold');
% legend('Location', 'north', 'FontSize', 10);

%% ========================================================================
%  4. THRUSTER SIZING – SLEW MANEUVER
% =========================================================================

fprintf('===========================================================\n');
fprintf('  4. THRUSTER SIZING – SLEW MANEUVER\n');
fprintf('===========================================================\n\n');

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

fprintf(' Acceleration/braking time  t_acc   = %.4f  s\n', t_acc);
fprintf(' Angle per phase            θ_acc   = %.6f  rad  (%.4f°)\n', theta_acc, theta_acc*(180/pi));
fprintf(' Coasting angle             θ_coast = %.4f  rad  (%.4f°)\n', theta_coast, theta_coast*(180/pi));
fprintf(' Coasting time              t_coast = %.2f  s\n', t_coast);
fprintf(' Total slew time            t_slew  = %.2f  s\n\n', t_slew_thr);

% -------------------------------------------------------------------------
%  Plot: Thruster slew angular rate profile
% -------------------------------------------------------------------------
t_thr_full = [0, t_acc, t_acc+t_coast, t_slew_thr];
r_thr      = [0, rate_max, rate_max, 0] * (180/pi);

% figure('Name','Thruster Slew Profile','Color','w','Position',[100 100 800 400]);
% hold on; grid on;
% plot(t_thr_full, r_thr, '-g', 'LineWidth', 2.5);
% yline(rate_max*(180/pi), '--k', 'LineWidth', 1.2, 'DisplayName','Rate limit 0.5°/s');
% xlabel('Time [s]', 'FontSize', 12);
% ylabel('Angular Rate [°/s]', 'FontSize', 12);
% title('Slew Maneuver – Angular Rate Profile (Thrusters)', 'FontSize', 13, 'FontWeight', 'bold');
% legend({'Thruster slew profile','Rate limit 0.5°/s'}, 'Location','northeast','FontSize',10);

%% ========================================================================
%  5. PROPELLANT MASS SIZING – WHEEL DESATURATION
% =========================================================================

fprintf('===========================================================\n');
fprintf('  5. PROPELLANT MASS SIZING – WHEEL DESATURATION\n');
fprintf('===========================================================\n\n');

% Theoretical minimum desaturation time (momentum transfer at full thrust)
t_desat_min = h_rw_max / (n_thr * L_thr * F_thr);
fprintf(' Theoretical minimum desaturation time = %.2f  s\n', t_desat_min);
fprintf(' (Requires RW torque >> T_max → physically unfeasible)\n\n');

% Realistic desaturation: driven by RW max torque (use T_max/2 conservatively)
% Time to dump h_max at T_max/2 from each side (momentum transferred by thrusters)
t_desat = 2 * h_rw_max / T_rw_max;

% Required force per thruster (couple: n = 2)
F_desat = h_rw_max / (n_thr * L_thr * t_desat);

fprintf(' Realistic desaturation (limited by T_max/2):\n');
fprintf(' t_desat  = 2·h_max / T_max = %.1f  s\n', t_desat);
fprintf(' F_desat  = h_max / (n·L·t) = %.4f  N\n\n', F_desat);

% Number of desaturation manoeuvres over 2-year mission
T_mission = 2 * 365.25 * 24 * 3600;   % [s]
N_desat   = (T_mission / T_orb) / N_orb;
fprintf(' Mission duration    = 2 years = %.2e  s\n', T_mission);
fprintf(' Total orbits        = %.0f\n', T_mission / T_orb);
fprintf(' Desaturation events = %.0f\n\n', N_desat);

% Propellant mass (Tsiolkovsky / rocket equation form for impulsive burns)
% M_prop = I_total / (Isp · g0),   I_total = N_desat · n · t_desat · F_desat
I_total  = N_desat * n_thr * t_desat * F_desat;
M_prop   = I_total / (Isp * g0);

fprintf(' Total impulse for desaturation:\n');
fprintf(' I_total = N_desat × n × t_desat × F_desat = %.2f  Ns\n', I_total);
fprintf(' Propellant mass M_prop = I_total / (Isp·g0) = %.2f  kg\n\n', M_prop);
fprintf(' !! Propulsion subsystem margins shall be applied on top of this value.\n\n');

%% ========================================================================
%  6. SUMMARY TABLE
% =========================================================================

fprintf('===========================================================\n');
fprintf('  6. RESULTS SUMMARY\n');
fprintf('===========================================================\n\n');

fprintf(' %-23s %13s %12.4e %6s\n', 'GG Torque',                "T_gg", T_gg, 'Nm');
fprintf(' %-23s %13s %12.4e %6s\n', 'SRP Torque',               "T_srp", T_srp, "Nm");
fprintf(' %-23s %13s %12.4e %6s\n', 'Tot. disturbance',         "T_dis", T_dis, "Nm");
fprintf(' %-23s %13s %12.4f %6s\n', 'Momentum per orbit',       "h_orbit", h_orbit, "Nms");
fprintf(' %-23s %13s %12.2f %6s\n', 'Orbits before RW satur.',  "N_orb", N_orb, "orb");
fprintf('\n');
fprintf(' %-23s %13s %12.2f %6s\n', 'Min slew time',            "t_slew_min", t_slew_min, "s");
fprintf(' %-23s %13s %12.4f %6s\n', 'Peak rate',                "-", rate_slew_max_A*(180/pi), "deg/s");
fprintf(' %-23s %13s %12.4f %6s\n', 'Reduced torque',           "T_reduced", T_reduced, "Nm");
fprintf(' %-23s %13s %12.1f %6s\n', 'Slew time',                "t_slew_red", t_slew_red, "s");
fprintf(' %-23s %13s %12.2f %6s\n', 'Peak slew RW momentum',    "h_slew_max_B", h_slew_max_B, "Nms");
fprintf('\n');
fprintf(' %-23s %13s %12.4f %6s\n', 'Thruster slew',            "t_acc", t_acc, "s");
fprintf(' %-23s %13s %12.2f %6s\n', 'Tot. Thruster slew time',  "t_slew_thr", t_slew_thr, "s");
fprintf('\n');
fprintf(' %-23s %13s %12.1f %6s\n', 'Real desat. duration',     "t_desat", t_desat, "s");
fprintf(' %-23s %13s %12.4f %6s\n', 'Thruster desat. force',    "F_desat", F_desat, "N");
fprintf(' %-23s %13s %12.0f %6s\n', 'Tot. desat. events 2 yr',  "N_desat", N_desat, "-");
fprintf(' %-23s %13s %12.2f %6s\n', 'Min propellant mass',      "M_prop", M_prop, "kg");
fprintf('\n===========================================================\n');
fprintf('  END OF ADCS SIZING\n');
fprintf('===========================================================\n');
