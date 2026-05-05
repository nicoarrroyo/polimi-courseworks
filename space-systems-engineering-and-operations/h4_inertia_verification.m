
clear; close all; clc;
fprintf('===========================================================\n');
fprintf('  ADCS SIZING  -  SMOS\n');
fprintf('===========================================================\n\n');

%% 1. KNOWN PARAMETERS (Masses in kg, Dimensions in meters)
% Proteus Bus (rectangular prism)
m_bus           = 225.08;
bus_lx          = 0.954;
bus_ly          = 0.954;
bus_lz          = 1.004;

% Solar Panel Arrays
% SMOS has two arrays, one on each side of the bus (+Y and -Y).
% Each array is a thin rectangular plate (lx ~ 0, ly = half-span, lz = height).
m_solar_total   = 49.92;                % Total mass of BOTH arrays [kg]
m_solar_one     = m_solar_total / 2;    % Mass of one array [kg]
solar_lx        = 0;                    % Panel thickness (assume zero) [m]
solar_ly        = 4.0;                  % Array span (along Y) [m]
solar_lz        = 1.5;                  % Panel height (along Z) [m]

% Y-offset of each array CoM from the satellite Z-axis:
% Assume massless strut connecting solar panels to bus surface with width 
% of 0.5 m todo-check. The solar array extends by half the bus width added 
% with this massless strut and half the length of the full array.
dy_solar    = (bus_ly / 2) + 0.5 + (solar_ly / 2);

% MIRAS Hub (solid cylinder, symmetry axis = Z)
m_hub           = 83.5;
hub_diam        = 1.3;          % Diameter of hub (hexagonal) [m]
hub_r           = hub_diam / 2;
hub_lz          = 1.2;          % Height of hub along Z [m]

% MIRAS Arms (3 arms, 120 deg apart, each assumed as a thin rod)
m_arm           = 90.5;         % Mass of a single arm [kg]
m_arms_total    = m_arm * 3;    % Total arms mass [kg]
arm_L           = 4.6;          % Length of each arm [m]
r_arm_start     = hub_r;        % Arm attachment radius (at hub surface) [m]

%% 2. Z-COORDINATES OF COMPONENT CoMs
% Datum: bottom face of hub = Z = 0.
% Stack order (bottom to top): Hub -> Solar -> Bus

z_arms_com  = 0; % Arm thickness is 0 todo-check
z_hub_com   = hub_lz / 2;
z_solar_com = hub_lz; % Solar panels attach at hub-bus interface todo-check
z_bus_com   = hub_lz + (bus_lz / 2);

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
    (m_bus/12) * (bus_ly^2 + bus_lz^2), ...  % Ixx
    (m_bus/12) * (bus_lx^2 + bus_lz^2), ...  % Iyy
    (m_bus/12) * (bus_lx^2 + bus_ly^2)  ...  % Izz
]);

% --- Hub (solid cylinder) ---
I_hub_local = diag([ ...
    (m_hub/12) * (3*hub_r^2 + hub_lz^2), ... % Ixx
    (m_hub/12) * (3*hub_r^2 + hub_lz^2), ... % Iyy
    (m_hub/2)  *    hub_r^2              ... % Izz
]);

% --- Solar Arrays (thin rectangular plates) ---
% Treat each array as an independent thin rectangular plate.
% Local inertia of ONE array about its own CoM (plate lies in YZ-plane, 
% thickness along X is negligible):
I_one_array_local = diag([ ...
    (m_solar_one/12) * (solar_ly^2 + solar_lz^2), ... % Ixx (width + height)
    (m_solar_one/12) * (solar_lx^2 + solar_lz^2), ... % Iyy (thickness + height)
    (m_solar_one/12) * (solar_lx^2 + solar_ly^2)  ... % Izz (thickness + width)
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
Izz_single_arm = (m_arm / 3) * (arm_L^2 + 3*r0*arm_L + 3*r0^2);

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
fprintf('   %-12s  Ixx=%10.2f  Iyy=%10.2f  Izz=%10.2f\n', 'Hub',   I_hub_shifted(1,1),   I_hub_shifted(2,2),   I_hub_shifted(3,3));
fprintf('   %-12s  Ixx=%10.2f  Iyy=%10.2f  Izz=%10.2f\n', 'Arms',  I_arms_shifted(1,1),  I_arms_shifted(2,2),  I_arms_shifted(3,3));
fprintf('   %-12s  Ixx=%10.2f  Iyy=%10.2f  Izz=%10.2f\n', 'Solar', I_solar_shifted(1,1), I_solar_shifted(2,2), I_solar_shifted(3,3));
fprintf('   %-12s  Ixx=%10.2f  Iyy=%10.2f  Izz=%10.2f\n', 'Bus',   I_bus_shifted(1,1),   I_bus_shifted(2,2),   I_bus_shifted(3,3));

fprintf('\n Total Inertia Tensor at Global CoM (kg*m^2):\n');
fprintf('   Ixx: %10.2f\n', I_total(1,1));
fprintf('   Iyy: %10.2f\n', I_total(2,2));
fprintf('   Izz: %10.2f\n', I_total(3,3));
fprintf('\n Note: Products of inertia are zero by 3-fold arm symmetry\n');
fprintf(' and array +/-Y symmetry about the XZ and YZ planes.\n');
fprintf('========================================\n\n');
