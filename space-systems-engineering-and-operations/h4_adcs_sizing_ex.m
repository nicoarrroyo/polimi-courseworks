
clear; close all; clc;
fprintf('===========================================================\n');
fprintf('  ADCS SIZING  -  SMOS\n');
fprintf('===========================================================\n\n');

%% 1. KNOWN PARAMETERS (Masses in kg, Dimensions in meters)
% Proteus Bus (rectangular prism)
m_bus           = 225.08;
lx_bus          = 0.954;
ly_bus          = 0.954;
lz_bus          = 1.004;

% Solar Panel Arrays
% SMOS has two arrays, one on each side of the bus (+Y and -Y).
% Each array is a thin rectangular plate (lx ~ 0, ly = half-span, lz = height).
m_solar_total   = 49.92;                % Total mass of BOTH arrays [kg]
m_solar_one     = m_solar_total / 2;    % Mass of one array [kg]
lx_solar        = 0;                    % Panel thickness (assume zero) [m]
ly_solar        = 4.0;                  % Array span (along Y) [m]
lz_solar        = 1.5;                  % Panel height (along Z) [m]

% Y-offset of each array CoM from the satellite Z-axis:
% Assume massless strut connecting solar panels to bus surface with width 
% of 0.5 m todo-check. The solar array extends by half the bus width added 
% with this massless strut and half the length of the full array.
dy_solar        = (ly_bus / 2) + 0.5 + (ly_solar / 2);

% MIRAS Hub (solid cylinder, symmetry axis = Z)
m_hub           = 83.5;
l_diam_hub        = 1.3;            % Diameter of hub (hexagonal) [m]
lr_hub           = l_diam_hub / 2;
lz_hub          = 1.2;              % Height of hub along Z [m]

% MIRAS Arms (3 arms, 120 deg apart, each assumed as a thin rod)
m_arm           = 90.5;         % Mass of a single arm [kg]
m_arms_total    = m_arm * 3;    % Total arms mass [kg]
l_arm           = 4.6;          % Length of each arm [m]
r_arm_start     = lr_hub;       % Arm attachment radius (at hub surface) [m]

% Total height
lz_smos = lz_bus + lz_hub;

%% 2. Z-COORDINATES OF COMPONENT CoMs
% Datum: bottom face of hub = Z = 0.
% Stack order (bottom to top): Hub -> Solar -> Bus

z_arms_com  = 0; % Arm thickness is 0 todo-check
z_hub_com   = lz_hub / 2;
z_solar_com = lz_hub; % Solar panels attach at hub-bus interface todo-check
z_bus_com   = lz_hub + (lz_bus / 2);

%% 3. GLOBAL CENTER OF MASS (Z-axis)
total_mass = m_bus + m_solar_total + m_hub + m_arms_total;

z_cg = ((m_bus          * z_bus_com)   + ...
        (m_solar_total  * z_solar_com) + ...
        (m_hub          * z_hub_com)   + ...
        (m_arms_total   * z_arms_com) ) ... 
        / total_mass;

%% 4. LOCAL INERTIA TENSORS (About each component's own CoM)
% --- Bus (rectangular prism) ---
I_bus_local = diag([ ...
    (m_bus/12) * (ly_bus^2 + lz_bus^2), ...  % Ixx
    (m_bus/12) * (lx_bus^2 + lz_bus^2), ...  % Iyy
    (m_bus/12) * (lx_bus^2 + ly_bus^2)  ...  % Izz
]);

% --- Hub (solid cylinder) ---
I_hub_local = diag([ ...
    (m_hub/12) * (3*lr_hub^2 + lz_hub^2), ... % Ixx
    (m_hub/12) * (3*lr_hub^2 + lz_hub^2), ... % Iyy
    (m_hub/2)  *    lr_hub^2              ... % Izz
]);

% --- Solar Arrays (thin rectangular plates) ---
% Treat each array as an independent thin rectangular plate.
% Local inertia of ONE array about its own CoM (plate lies in YZ-plane, 
% thickness along X is negligible):
I_one_array_local = diag([ ...
    (m_solar_one/12) * (ly_solar^2 + lz_solar^2), ... % Ixx (width + height)
    (m_solar_one/12) * (lx_solar^2 + lz_solar^2), ... % Iyy (thickness + height)
    (m_solar_one/12) * (lx_solar^2 + ly_solar^2)  ... % Izz (thickness + width)
]);

% --- MIRAS Arms (thin rods todo-check-assumption) --- %
% Each arm is a thin rod lying in the XY-plane (0 thickness), starting at
% r0 = hub_r from the Z-axis and extending radially to r0 + arm_L.
%
% The three arms are at phi = 0, 120, 240 deg.
%
% Inertia of ONE arm, computed about the SATELLITE Z-axis directly:
%   Izz_arm  = integral_{r0}^{r0+L} r^2 * (m/L) dr
%             = (m/3L) * [(r0+L)^3 - r0^3]
%             = (m/3)  * (L^2 + 3*r0*L + 3*r0^2)
r0 = r_arm_start;
Izz_single_arm = (m_arm / 3) * (l_arm^2 + 3*r0*l_arm + 3*r0^2);

% Ixx and Iyy of one arm about satellite Z-axis:
% For a rod at angle phi in the XY-plane, each mass element dm at radius r
% has coordinates (r*cos(phi), r*sin(phi), 0).
%   Ixx_arm(phi) = integral r^2 * sin^2(phi) dm  =  sin^2(phi) * Izz_single_arm
%   Iyy_arm(phi) = integral r^2 * cos^2(phi) dm  =  cos^2(phi) * Izz_single_arm
% Summing over phi = 0, 120, 240 deg:
%   sum(sin^2) = 0 + 3/4 + 3/4 = 3/2
%   sum(cos^2) = 1 + 1/4 + 1/4 = 3/2
% Therefore by symmetry: Ixx_total = Iyy_total = (3/2) * Izz_single_arm
% NOTE: this numerical coincidence (Ixx=Iyy=1.5*Izz_single) is the same
% value as the original code, but now derived correctly.
%
% Izz_total = sum of all three Izz_arm(phi):
%   Each arm's Izz about satellite axis = Izz_single (independent of phi).
%   Total: 3 * Izz_single_arm
%
% Z-axis inertia: arms lie flat in XY-plane, so Ixx/Iyy about the arm CoM
% have no Z-extent contribution. The approach above already computes
% inertia about the SATELLITE Z-axis directly (no local+PAT split needed
% for the radial direction). We only need a Z-shift for the arm CoM
% sitting at z_arms_com vs z_cg.

dz_arms = z_arms_com - z_cg;

% Arms inertia about satellite Z-axis
I_arms_about_satZ = diag([ ...
    1.5 * Izz_single_arm, ...  % Ixx: sum sin^2(phi_k) * Izz_single
    1.5 * Izz_single_arm, ...  % Iyy: sum cos^2(phi_k) * Izz_single
    3.0 * Izz_single_arm  ...  % Izz: 3 * Izz_single
]);

%% 5. SHIFT COMPONENTS WITH PARALLEL AXIS THEOREM
% General Z-only shift (for axially symmetric components with CoM on Z-axis)
shift_z = @(mass, dz) diag([mass*dz^2, mass*dz^2, 0]);

% --- Bus ---
dz_bus = z_bus_com - z_cg;
I_bus_shifted = I_bus_local + shift_z(m_bus, dz_bus);

% --- Hub ---
dz_hub = z_hub_com - z_cg;
I_hub_shifted = I_hub_local + shift_z(m_hub, dz_hub);

% --- Solar Arrays ---
% Parallel axis theorem for each array: shift by (dy_solar, 0, dz_solar).
% array +Y: offset (+dy_solar, 0, dz_solar_from_cg)
% array -Y: offset (-dy_solar, 0, dz_solar_from_cg)  [dy^2 is same]
dz_solar = z_solar_com - z_cg;
% Shift tensor for a general (dx, dy, dz) offset:
%   diag([dy^2+dz^2, dx^2+dz^2, dx^2+dy^2]) * mass
% For a Y-offset array, dx = 0:
shift_one_array = @(dy, dz) diag([ ...
    m_solar_one * (dy^2 + dz^2), ... % Ixx
    m_solar_one * (0   + dz^2), ...  % Iyy  (no X offset)
    m_solar_one * (dy^2 + 0  )  ...  % Izz  (no Z contribution to Izz)
]);
% Both arrays have the same |dy| and same dz, so sum = 2 * shifted_one_array
I_solar_shifted = 2 * (I_one_array_local + shift_one_array(dy_solar, dz_solar));

% --- MIRAS Arms ---
% Apply Z-shift from arm plane (z_arms_com) to global CoM:
% This adds m_arms_total * dz^2 to Ixx and Iyy only.
I_arms_shifted = I_arms_about_satZ + diag([ ...
    m_arms_total * dz_arms^2, ...
    m_arms_total * dz_arms^2, ...
    0 ...
]);

%% 6. TOTAL INERTIA TENSOR AT GLOBAL CoM
I_total = I_bus_shifted + I_solar_shifted + I_hub_shifted + I_arms_shifted;

%% 7. DISPLAY RESULTS

fprintf('========================================\n');
fprintf(' SMOS SATELLITE INERTIA\n');
fprintf('========================================\n');
fprintf(' Total System Mass         : %.2f kg\n', total_mass);
fprintf(' Global CoM (Z from datum) : %.4f m\n\n', z_cg);

fprintf(' Component Z-offsets from global CoM:\n');
fprintf('   Hub    dz = %+.4f m\n', dz_hub);
fprintf('   Arms   dz = %+.4f m\n', dz_arms);
fprintf('   Solar  dz = %+.4f m  (lateral offset dy = %.4f m)\n', dz_solar, dy_solar);
fprintf('   Bus    dz = %+.4f m\n\n', dz_bus);

fprintf(' Inertia contributions (kg*m^2):\n');
fprintf('   %-6s  Ixx=%10.2f | Iyy=%10.2f | Izz=%10.2f\n', 'Hub',   ...
    I_hub_shifted(1,1),   I_hub_shifted(2,2),   I_hub_shifted(3,3));
fprintf('   %-6s  Ixx=%10.2f | Iyy=%10.2f | Izz=%10.2f\n', 'Arms',  ...
    I_arms_shifted(1,1),  I_arms_shifted(2,2),  I_arms_shifted(3,3));
fprintf('   %-6s  Ixx=%10.2f | Iyy=%10.2f | Izz=%10.2f\n', 'Solar', ...
    I_solar_shifted(1,1), I_solar_shifted(2,2), I_solar_shifted(3,3));
fprintf('   %-6s  Ixx=%10.2f | Iyy=%10.2f | Izz=%10.2f\n', 'Bus',   ...
    I_bus_shifted(1,1),   I_bus_shifted(2,2),   I_bus_shifted(3,3));

fprintf('\n Total Inertia Tensor at Global CoM (kg*m^2):\n');
fprintf('   Ixx: %10.2f\n', I_total(1,1));
fprintf('   Iyy: %10.2f\n', I_total(2,2));
fprintf('   Izz: %10.2f\n', I_total(3,3));
fprintf('\n Note: Products of inertia are zero by 3-fold arm symmetry\n');
fprintf(' and array +/-Y symmetry about the XZ and YZ planes.\n');
fprintf('========================================\n\n');

%% ========================================================================
%  1. MISSION & ENVIRONMENT PARAMETERS
% =========================================================================

% --- Spacecraft parameters -----------------------------------------------
Ix          = I_total(1,1);         % Moment of inertia – X axis [kg m²]
Iy          = I_total(2,2);         % Moment of inertia – Y axis [kg m²]
Iz          = I_total(3,3);         % Moment of inertia – Z axis [kg m²]
I_all       = [Ix, Iy, Iz];
I_min       = min(I_all);           % Minimum inertia [kg m²]
I_max       = max(I_all);           % Maximum inertia [kg m²]

S_solar     = 2*ly_solar*lz_solar;  % Solar array cross-sectional area
S_bus       = lz_bus * ly_bus;      % PROTEUS bus cross-sectional area
S_hub       = lz_hub * l_diam_hub;  % MIRAS hub cross-sectional area
S           = S_solar+S_bus+S_hub;  % Cross-sectional area [m²]
d_com       = (lz_smos / 2) - z_cg; % CoM displacement (csp – cg) [m]
q_refl      = 0.3;                  % Reflectivity coefficient [-] TODO

T_rw_max    = 0.075;                % RW max torque [Nm] TODO
h_rw_max    = 12.0;                 % RW max momentum [Nms] TODO

F_thr       = 10;                   % Thruster force [N] TODO
L_thr       = 0.69;                 % Thruster moment arm [m] TODO
n_thr       = 4;                    % Number of thrusters in couple [-] TODO
Isp         = 292;                  % Specific impulse [s] TODO
g0          = 9.81;                 % Standard gravity [m/s²]

rate_max    = 0.5 * (pi/180);       % Max angular rate [rad/s] TODO
theta_max   = 180 * (pi/180);       % Worst-case slew angle [rad] TODO

% --- Orbital / environment parameters ------------------------------------
mu      = 398600.4e9;   % Earth gravitational constant  [m³/s²]
R_E     = 6371e3;       % Earth equatorial radius       [m]
a       = 765e3 + R_E;  % Semi-major axis               [m]
T_orb   = 2*pi * sqrt((a^3) / mu); % Orbit period                  [s]
F_sun   = 1362;         % Solar constant at Earth       [W/m²]
c_light = 3e8;          % Speed of light                [m/s]

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

T_srp = (F_sun / c_light) * S * (1 + q_refl) * cos(i_sun) * d_com;

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
