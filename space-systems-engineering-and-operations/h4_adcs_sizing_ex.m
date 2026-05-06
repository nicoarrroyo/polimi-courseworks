clear; close all; clc;
fprintf('===========================================================\n');
fprintf('  ADCS SIZING  -  SMOS\n');
fprintf('===========================================================\n\n');

%% 1. KNOWN PARAMETERS (Masses in kg, Dimensions in meters)
% Proteus Bus (rectangular prism)
m_bus           = 225.08+28;
lx_bus          = 1.004;
ly_bus          = 0.954;
lz_bus          = 0.954;

% Solar Panel Arrays
% SMOS has two arrays, one on each side of the bus (+Y and -Y).
% Each array is a thin rectangular plate (lz ~ 0, ly = length, lx = height).
m_solar_total   = 49.92;                % Total mass of BOTH arrays [kg]
m_solar         = m_solar_total / 2;    % Mass of one array [kg]
lx_solar        = 1.5;                  % Panel height [m]
ly_solar        = 4.0;                  % Array span (along Y) [m]
lz_solar        = 0;                    % Panel Thickness (negligible) [m]

% Y-offset of each array CoM from the satellite Z-axis:
% Assume massless strut connecting solar panels to bus surface with width 
% of 0.5 m todo-check. The solar array extends by half the bus width added 
% with this massless strut and half the length of the full array.
dy_solar        = (ly_bus + ly_solar)/2 + 0.5;

% MIRAS Arms (3 arms, 120 deg apart, each assumed as a thin rod)
m_arm           = 64;           % Mass of a single arm [kg]
m_arms_total    = m_arm * 3;    % Total arms mass [kg]
L_arm           = 3.4;          % Length of each arm [m]
r_arm_start     = 1.3/2;        % Arm attachment radius (at hub surface) [m]

% MIRAS Hub (solid cylinder, symmetry axis = X)
m_hub           = 355 - m_arms_total;
hub_diam        = 1.3;          % Diameter of hub (hexagonal) [m]
lr_hub          = hub_diam / 2;
lx_hub          = 1.2;          % Height of hub along X [m]

%% 2. Z-COORDINATES OF COMPONENT CoMs
% Datum: bottom face of hub = Z = 0.
% Stack order (bottom to top): Hub -> Solar -> Bus

x_arms_com  = lx_hub + lx_bus;          % Arm thickness negligible 
x_hub_com   = lx_bus + (lx_hub / 2);
x_solar_com = lx_bus*0.9;               % Solar panels attach near hub-bus interface
x_bus_com   = lx_bus/2;

%% 3. GLOBAL CENTER OF MASS (X-axis)
total_mass = m_bus + m_solar_total + m_hub + m_arms_total;

x_cg = ((m_bus          * x_bus_com)   + ...
        (m_solar_total  * x_solar_com) + ...
        (m_hub          * x_hub_com)   + ...
        (m_arms_total   * x_arms_com) ) ... 
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
    (m_hub/2)  *  lr_hub^2,               ... % Ixx
    (m_hub/12) * (3*lr_hub^2 + lx_hub^2), ... % Iyy
    (m_hub/12) * (3*lr_hub^2 + lx_hub^2), ... % Izz
]);

% --- Solar Arrays (thin rectangular plates) ---
% Treat each array as an independent thin rectangular plate.
% Local inertia of ONE array about its own CoM (plate lies in YZ-plane, 
% thickness along X is negligible):
I_array_local = diag([ ...
    (m_solar/12) * (ly_solar^2 + lz_solar^2), ... % Ixx (width + height)
    (m_solar/12) * (lx_solar^2 + lz_solar^2), ... % Iyy (thickness + height)
    (m_solar/12) * (lx_solar^2 + ly_solar^2)  ... % Izz (thickness + width)
]);

% --- MIRAS Arms (thin rods todo-check-assumption) --- %
% Each arm is a thin rod lying in the XY-plane (0 thickness), starting at
% r0 = hub_r from the Z-axis and extending radially to r0 + arm_L.
%
% The three arms are at phi = 0, 120, 240 deg.

Ixx_single_arm = (m_arm / 3)*L_arm^2 + m_arm*r_arm_start^2;

% Ixx and Iyy of one arm about satellite Z-axis:
% For a rod at angle phi in the XY-plane, each mass element dm at radius r
% has coordinates (r*cos(phi), r*sin(phi), 0).
%   Izz_arm(phi) = integral r^2 * sin^2(phi) dm  =  sin^2(phi) * Ixx_single_arm
%   Iyy_arm(phi) = integral r^2 * cos^2(phi) dm  =  cos^2(phi) * Ixx_single_arm
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

dx_arms = x_arms_com - x_cg;

% Arms inertia tensor
I_arms = diag([ ...
    3.0 * Ixx_single_arm, ...  % Ixx: 3 * Ixx_single
    1.5 * Ixx_single_arm, ...  % Iyy: sum cos^2(phi_k) * Ixx_single
    1.5 * Ixx_single_arm  ...  % Izz: sum sin^2(phi_k) * Ixx_single
]);

%% 5. SHIFT COMPONENTS WITH PARALLEL AXIS THEOREM
shift_x = @(mass, dx) mass*diag([0, dx^2, dx^2]);
shift_y = @(mass, dy) mass*diag([dy^2, 0, dy^2]);

% --- Bus ---
dx_bus = x_bus_com - x_cg;
I_bus_shifted = I_bus_local + shift_x(m_bus, dx_bus);

% --- Hub ---
dx_hub = x_hub_com - x_cg;
I_hub_shifted = I_hub_local + shift_x(m_hub, dx_hub);

% --- Solar Arrays ---
dx_solar = x_solar_com - x_cg;

% Both arrays have the same |dy| and same dx, so sum = 2 * shifted_one_array
I_solar_shifted = 2*(I_array_local + shift_x(m_solar, dx_solar) + shift_y(m_solar, dy_solar));

% --- MIRAS Arms ---
% Apply Z-shift from arm plane (z_arms_com) to global CoM:
% This adds m_arms_total * dz^2 to Ixx and Iyy only.
I_arms_shifted = I_arms + shift_x(m_arms_total, dx_arms);

%% 6. TOTAL INERTIA TENSOR AT GLOBAL CoM
I_total = I_bus_shifted + I_solar_shifted + I_hub_shifted + I_arms_shifted;

%% 7. DISPLAY RESULTS

fprintf('========================================\n');
fprintf(' SMOS SATELLITE INERTIA\n');
fprintf('========================================\n');
fprintf(' Total System Mass         : %.2f kg\n', total_mass);
fprintf(' Global CoM (X from datum) : %.4f m\n\n', x_cg);

fprintf(' Component Z-offsets from global CoM:\n');
fprintf('   Arms   dx = %+.4f m\n', dx_arms);
fprintf('   Hub    dx = %+.4f m\n', dx_hub);
fprintf('   Solar  dx = %+.4f m  (lateral offset dy = %.4f m)\n', dx_solar, dy_solar);
fprintf('   Bus    dx = %+.4f m\n\n', dx_bus);

fprintf(' Inertia contributions LOCAL (kg*m^2):\n');
fprintf('   %-6s  Ixx=%10.2f | Iyy=%10.2f | Izz=%10.2f\n', 'Arms',  ...
    I_arms(1,1),  I_arms(2,2),  I_arms(3,3));
fprintf('   %-6s  Ixx=%10.2f | Iyy=%10.2f | Izz=%10.2f\n', 'Hub',   ...
    I_hub_local(1,1),   I_hub_local(2,2),   I_hub_local(3,3));
fprintf('   %-6s  Ixx=%10.2f | Iyy=%10.2f | Izz=%10.2f\n', 'Solar', ...
    2*I_array_local(1,1), 2*I_array_local(2,2), 2*I_array_local(3,3));
fprintf('   %-6s  Ixx=%10.2f | Iyy=%10.2f | Izz=%10.2f\n', 'Bus',   ...
    I_bus_local(1,1),   I_bus_local(2,2),   I_bus_local(3,3));

fprintf(' Inertia contributions SHIFTED (kg*m^2):\n');
fprintf('   %-6s  Ixx=%10.2f | Iyy=%10.2f | Izz=%10.2f\n', 'Arms',  ...
    I_arms_shifted(1,1),  I_arms_shifted(2,2),  I_arms_shifted(3,3));
fprintf('   %-6s  Ixx=%10.2f | Iyy=%10.2f | Izz=%10.2f\n', 'Hub',   ...
    I_hub_shifted(1,1),   I_hub_shifted(2,2),   I_hub_shifted(3,3));
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

%% 7. CENTER OF PRESSURE for SRP

A_panel = lx_solar*ly_solar;
A_bus = lx_bus*ly_bus;
A_hub = lx_hub*hub_diam;

r_panel = x_solar_com;
r_bus = x_bus_com;
r_hub = x_hub_com;

Q_panel = 1.1;
Q_bus = 1.3;
Q_hub = 1.3;

x_cp = (2*A_panel*Q_panel*r_panel + A_bus*Q_bus*r_bus + A_hub*Q_hub*r_hub) / ...
    (2*A_panel*Q_panel+ A_bus*Q_bus + A_hub*Q_hub);

csp_shift = x_cg - x_cp;

%% ========================================================================
%  1. MISSION & ENVIRONMENT PARAMETERS
% =========================================================================

% --- Spacecraft parameters -----------------------------------------------
Ix      = I_total(1,1);
Iy      = I_total(2,2);
Iz      = I_total(3,3);
I       = [Ix, Iy, Iz]; % Full inertia tensor (as a row)
I_min   = min(I);       % Minimum inertia [kg m²]
I_max   = max(I);       % Maximum inertia [kg m²]

S       = 2*A_panel + A_bus + A_hub;        % Cross-sectional area [m²]
q_panel = 0.1*2*A_panel;
q_bus   = 0.3*A_bus;
q_hub   = 0.3*A_hub;
q_refl  = (q_panel + q_bus + q_hub) / S;    % Reflectivity coefficient [-]

T_rw_max    = 0.120;            % RW max torque         [Nm]
T_max = T_rw_max * sqrt(3);
h_rw_max    = 8.00;             % RW max momentum       [Nms]
h_max = h_rw_max * sqrt(3);

rate_max    = 0.1 * (pi/180);   % Max angular rate      [rad/s] TODO
theta_max   = 180 * (pi/180);   % Worst-case slew angle [rad] TODO

% --- Orbital / environment parameters ------------------------------------
mu      = 398600.4e9;           % Earth gravitational constant  [m³/s²]
R_E     = 6371e3;               % Earth equatorial radius       [m]
a       = 765e3 + R_E;          % Semi-major axis               [m]
T_orb   = 2*pi*sqrt((a^3)/mu);  % Orbit period                  [s]
F_sun   = 1362;                 % Solar constant at Earth       [W/m²]
c_light = 3e8;                  % Speed of light                [m/s]

% Orbital radius (approximated as semi-major axis – conservative)
R_orb     = a;

fprintf('--- Spacecraft Parameters ---\n');
fprintf(' Inertia:  Ix = %.2f kg m^2\n', Ix);
fprintf(' Inertia:  Iy = %.2f kg m^2\n', Iy);
fprintf(' Inertia:  Iz = %.2f kg m^2\n', Iz);
fprintf(' Inertia:  I_min = %.2f kg m^2\n', I_min);
fprintf(' Inertia:  I_max = %.2f kg m^2\n', I_max);
fprintf(' RW:       max torque = %.3f Nm\n', T_rw_max);
fprintf(' RW:       max momentum = %.1f Nms\n', h_rw_max);

%% ========================================================================
%  2. DISTURBANCE ENVIRONMENT EVALUATION
% =========================================================================

fprintf('===========================================================\n');
fprintf('  2. DISTURBANCE ENVIRONMENT\n');
fprintf('===========================================================\n\n');
% Done by Preliminary Disturbances Evaluation

% --- Gravity Gradient ---
T_gg = 3*mu / (2*R_orb^3) * abs(I_max - I_min);

% --- Solar Radiation Pressure ---)
T_srp = (F_sun/c_light) * S * (1+q_refl) * (csp_shift);

% --- Magnetic Field ---
B_N_max = 4e-5;
B_N_min = 2e-5;
D_res = 1;
T_mag = D_res * B_N_max;

% --- Aero drag ---
rho = 1.25 * 10^-14;
V = sqrt(mu / a);
S_aero_arms = L_arm * sind(32.5) * 0.4 * 3;
S_aero_hub_lower_face = pi * (lr_hub)^2 * sind(32.5);
S_aero_hub_body = lx_hub * hub_diam * cosd(32.5);
S_aero_bus_body = lx_bus * ly_bus * cosd(32.5);
A_drag = S_aero_arms + S_aero_hub_lower_face + S_aero_bus_body + S_aero_hub_body;
x_cp_aero = ( ...
    (S_aero_arms+S_aero_hub_lower_face)*x_arms_com + ...
    S_aero_bus_body*x_bus_com + ...
    S_aero_hub_body*x_hub_com) / ...
    (A_drag);
cp_shift = x_cg - x_cp_aero;
Cd = 2.6; % from fauste
T_aero = 0.5 * rho * V^2 * A_drag * Cd * cp_shift;

fprintf(' T_gg = %.4e Nm\n', T_gg);
fprintf(' T_srp = %.4e Nm\n', T_srp);
fprintf(' T_mag = %.4e Nm\n', T_mag);
fprintf(' T_aero = %.4e Nm\n\n', T_aero);

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

T_dis = 2 * (T_gg+T_srp+T_aero);  % Total constant disturbance torque [Nm]
h_orbit = T_dis * T_orb;          % Momentum accumulated per orbit [Nms]

% Number of orbits before RW saturation
N_orb = h_max / h_orbit;

fprintf('--- 3.1  Disturbance Momentum Accumulation ---\n');
fprintf(' T_dis = 2 × T_gg = %.4e Nm\n', T_dis);
fprintf(' h_orbit = T_dis × T_orb = %.4f Nms\n', h_orbit);
fprintf(' N_orb (before saturation) = h_max / h_orbit = %.2f orbits\n\n', N_orb);

% Check cyclic Mag peak (should be negligible)
omega_orb    = 2*pi / T_orb; % [rad s^-1]
h_cyclic_max = (T_mag / omega_orb) * 2;   % integral of sin over half-period × 2
fprintf(' Cyclic MAG peak momentum = %.4f  Nms  (≪ h_max = %.1f Nms ✓)\n\n', h_cyclic_max, h_rw_max);

% -------------------------------------------------------------------------
%  Plot: Momentum accumulation over n orbits
% -------------------------------------------------------------------------
t_plot   = linspace(0, N_orb * T_orb, 5000);
h_sec    = T_dis * t_plot; % Secular (GG)

% -------------------------------------------------------------------------
%  3.2  Slew Maneuver Sizing
% -------------------------------------------------------------------------
fprintf('--- 3.2  Slew Maneuver ---\n');
fprintf(' Worst-case slew angle : θ_slew = %.0f°\n', theta_max*(180/pi));

% Minimum slew time using maximum torque
t_slew_min = sqrt(4 * theta_max * I_max / T_max);
fprintf('\n  [Case A] Maximum torque applied:\n');
fprintf(' t_slew_min = sqrt(4·θ·I_max / T_max) = %.2f  s\n', t_slew_min);

% Peak angular rate and momentum at max torque
rate_slew_max_A = (T_max / I_max) * (t_slew_min / 2);
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

%% ========================================================================
%  4. MAGNETORQUER SIZING – SLEW MANEUVER
% =========================================================================

fprintf('===========================================================\n');
fprintf('  4. MAGNETORQUER SIZING – SLEW MANEUVER\n');
fprintf('===========================================================\n\n');

D_torquer_min = T_dis / B_N_min;

%% ========================================================================
%  5. SUMMARY TABLE
% =========================================================================

fprintf('===========================================================\n');
fprintf('  5. RESULTS SUMMARY\n');
fprintf('===========================================================\n\n');

fprintf(' %-23s %13s %12.4e %6s\n', 'GG Torque',            "T_gg", T_gg, 'Nm');
fprintf(' %-23s %13s %12.4e %6s\n', 'SRP Torque',           "T_srp", T_srp, "Nm");
fprintf(' %-23s %13s %12.4e %6s\n', 'Mag Torque',           "T_mag", T_mag, "Nm");
fprintf(' %-23s %13s %12.4e %6s\n', 'Aero Torque',          "T_aero", T_aero, "Nm");
fprintf(' %-23s %13s %12.4e %6s\n', 'Tot. disturbance',     "T_dis", T_dis, "Nm");
fprintf(' %-23s %13s %12.4f %6s\n', 'Momentum per orbit',   "h_orbit", h_orbit, "Nms");
fprintf(' %-23s %13s %12.2f %6s\n', 'RW satur. in',         "N_orb", N_orb, "orbits");
fprintf('\n');
fprintf(' %-23s %13s %12.2f %6s\n', 'Min slew time',        "t_slew_min", t_slew_min, "s");
fprintf(' %-23s %13s %12.4f %6s\n', 'Peak rate',            "-", rate_slew_max_A*(180/pi), "deg/s");
fprintf(' %-23s %13s %12.4f %6s\n', 'Reduced torque',       "T_reduced", T_reduced, "Nm");
fprintf(' %-23s %13s %12.1f %6s\n', 'Slew time',            "t_slew_red", t_slew_red, "s");
fprintf(' %-23s %13s %12.2f %6s\n', 'Peak slew RW momentum',"h_slew_max_B", h_slew_max_B, "Nms");
fprintf('\n');
fprintf('\n===========================================================\n');
fprintf('  END OF ADCS SIZING\n');
fprintf('===========================================================\n');
