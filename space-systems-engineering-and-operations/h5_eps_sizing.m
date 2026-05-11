%% EPS Sizing
% This script sizes the Electrical Power System (EPS) of a spacecraft in
% Sun-Synchronous Orbit (SSO), covering:
%   1. Solar array sizing (area, cell count, series/parallel configuration)
%   2. Battery sizing (energy capacity, cell count, series/parallel config)

clc; clear; close all;

%% =========================================================================
%  SECTION 0: MISSION & ORBIT PARAMETERS
%  =========================================================================

R_E     = 6371;                 % Earth radius [km]
h       = 765;                  % Altitude [km]
mu      = 398600.4;             % Earth gravitational parameter [km^2 s^{-3}]
a       = h + R_E;              % Semi-major axis
e       = 0.0001;               % Eccentricity [-] (negligible)
i       = 98.44;                % Inclination [deg]
T_orb   = 2*pi * sqrt(a^3/mu);  % Orbit period [sec]

T_orb_min = T_orb / 60;         % Orbital period [min]
T_orb_h   = T_orb_min / 60;     % Orbital period [h]

beta_deg = 0;                   % Sun-orbit plane angle [deg] — set worst case
beta_max_deg = asind(R_E / a);  % Max beta for eclipse to occur

if abs(beta_deg) >= beta_max_deg
    F_ecl = 0;
else
    arg = sqrt(1 - (R_E/a)^2) / cosd(beta_deg);
    F_ecl = acosd(arg) / 180;   % Eclipse fraction of orbit
end

T_ecl_min = F_ecl * T_orb_min; % Eclipse duration [h]
T_sun_min  = T_orb_min - T_ecl_h;

% Eclipsed orbit case (during winter eclipse period)
T_ecl_win   = 20*60;                % Time in eclipse per orbit [sec] TODO-calculate
T_sun_win   = T_orb - T_ecl_win;    % Time in sunlight per orbit [sec]

% Nominal orbit case
T_ecl_nom   = 0;
T_sun_nom   = T_orb;

fprintf('=== ORBIT PARAMETERS ===\n');
fprintf(' Semi-Major Axis:      %g km (SSO, LTAN 06:00)\n', a);
fprintf(' Orbital period:       %.2f seconds\n', T_orb);
fprintf(' Winter sunlit time:   %.2f minutes\n', T_sun_win/60);
fprintf(' Winter eclipse time:  %.2f minutes\n', T_ecl_win/60);
fprintf(' Nominal sunlit time:  %.2f minutes\n', T_sun_nom/60);
fprintf(' Nominal eclipse time: %.2f minutes\n\n', T_ecl_nom/60);

%% =========================================================================
%  SECTION 1: POWER BUDGET (from Table 2)
%  =========================================================================
% Modes:
% Safe-Hold, Star Acquisition, Normal Autonomous, Orbit Correction 2/4
% All values in [W], including 20% system margin

% Total required power in Safe-Hold Mode [W]
P_shm       = 13.82*1.2;

% Total required power in Star Acquisition Mode [W]
P_sam       = 144.82*1.2;

% Total required power in Normal Autonomous Mode nominally [W]
P_nam_nom   = 143.70*1.2;

% Total required power in Normal Autonomous Mode during eclipse [W]
P_nam_ecl   = P_nam_nom*1.5*1.2; % TODO

% Total required power in Orbit COrrection Mode 2 [W]
P_ocm2      = 162.70*1.2;

% Total required power in Orbit Correction Mode 4 [W]
P_ocm4      = 181.70*1.2;

fprintf('=== POWER BUDGET (incl. 20%% margin) ===\n');
fprintf(' Science-Sunlight: %g W\n', P_science_sun);
fprintf(' Eclipse:          %g W\n', P_eclipse);
fprintf(' Telecom:          %g W\n', P_telecom);
fprintf(' Safe Mode:        %g W\n\n', P_safe);

%% =========================================================================
%  SECTION 2: EPS CONFIGURATION & EFFICIENCIES TODO
%  =========================================================================

control_type = 'PPT';   % Power conditioning type: 'PPT' or 'DET'

% Line efficiencies (Table 3)
eff_PPT_ecl = 0.60;  eff_PPT_sun = 0.80;
eff_DET_ecl = 0.65;  eff_DET_sun = 0.85;

switch upper(control_type)
    case 'PPT'
        X_e = eff_PPT_ecl;
        X_s = eff_PPT_sun;
    case 'DET'
        X_e = eff_DET_ecl;
        X_s = eff_DET_sun;
    otherwise
        error('Unknown control type: %s. Use ''PPT'' or ''DET''.', control_type);
end

fprintf('=== EPS CONFIGURATION ===\n');
fprintf(' Control type:         %s\n', control_type);
fprintf(' Eclipse efficiency Xe: %.2f\n', X_e);
fprintf(' Sunlight efficiency Xs: %.2f\n\n', X_s);

%% =========================================================================
%  SECTION 3: SOLAR ARRAY SIZING
%  =========================================================================

fprintf('=== SOLAR ARRAY SIZING ===\n');

% --- 3.1  Required power from solar arrays at EoL (energy balance, Eq. 1) ---
%
%         (Pe * Te / Xe)  +  (Ps * Ts / Xs)
%  Psa = ----------------------------------------
%                        Ts

P_sa = ( (P_eclipse * T_ecl_h / X_e) + (P_science_sun * T_sun_min / X_s) ) ...
       / T_sun_min;

fprintf('  Required solar array power at EoL:\n');
fprintf('    Psa = (Pe*Te/Xe + Ps*Ts/Xs) / Ts\n');
fprintf('        = (%.0f*%.2f/%.2f + %.0f*%.2f/%.2f) / %.2f\n', ...
    P_eclipse, T_ecl_h, X_e, P_science_sun, T_sun_min, X_s, T_sun_min);
fprintf('        = %.1f W\n\n', P_sa);

% --- 3.2  Specific power at BoL (Eq. 2) ---
%
%  Pbol = P0 * epsilon * Id * cos(theta)

P0          = 1365;     % Solar constant at Earth [W m^{-2}]
epsilon     = 0.15;     % Silicon cell efficiency [-]
Id          = 0.80;     % Inherent degradation factor [-]
SAA_deg     = 5;        % Sun Aspect Angle [degrees]
theta_rad   = deg2rad(SAA_deg);

P_bol_spec = P0 * epsilon * Id * cos(theta_rad);  % [W m^{-2}]

fprintf('  Specific power at BoL (Eq. 2):\n');
fprintf('    P0 = %g W/m²,  epsilon = %.2f,  Id = %.2f,  SAA = %g deg\n', ...
    P0, epsilon, Id, SAA_deg);
fprintf('    Pbol = P0 * eps * Id * cos(SAA) = %.1f W/m²\n\n', P_bol_spec);

% --- 3.3  Specific power at EoL (Eq. 3) ---
%
%  Peol = Pbol * (1 - dpy)^mission_years

dpy          = 0.03;    % Degradation per year [-]
mission_yrs  = 12;      % Mission duration [years]

P_eol_spec = P_bol_spec * (1 - dpy)^mission_yrs;  % [W m^{-2}]

fprintf('  Specific power at EoL (Eq. 3):\n');
fprintf('    dpy = %.2f,  mission duration = %g years\n', dpy, mission_yrs);
fprintf('    Peol = Pbol * (1 - dpy)^n = %.1f W/m²\n\n', P_eol_spec);

% --- 3.4  Required solar array area (Eq. 4) ---
A_sa = P_sa / P_eol_spec;   % [m^2]

fprintf('  Required solar array area:\n');
fprintf('    Asa = Psa / Peol = %.1f / %.1f = %.2f m²\n\n', ...
    P_sa, P_eol_spec, A_sa);

% --- 3.5  Cell sizing (CESI CTJ30 reference cell) ---
A_cell   = 30.15e-4;    % Cell area [m^2]  (30.15 cm²)
V_cell   = 2.33;        % Cell max-power voltage [V]
V_bus_min    = 23;          % Bus voltage [V]
V_bus_max    = 36;          % Bus voltage [V]

% Minimum total number of cells
N_cells_min = ceil(A_sa / A_cell);

% Number of cells in series (to match bus voltage)
N_series = ceil(V_bus_min / V_cell);

% Number of parallel strings (+ 1 redundant string)
N_parallel = ceil(N_cells_min / N_series) + 1;

% Actual total cell count and array area
N_cells_real = N_parallel * N_series;
A_sa_real    = A_cell * N_cells_real;   % [m^2]

fprintf('  Cell configuration (CESI CTJ30: A=%.2f cm², Vcell=%.2f V):\n', ...
    A_cell*1e4, V_cell);
fprintf('    Minimum cells required:  %d\n', N_cells_min);
fprintf('    Cells in series:         %d  (Vbus/Vcell = %g/%.2f)\n', ...
    N_series, V_bus_min, V_cell);
fprintf('    Parallel strings:        %d  (incl. 1 redundant string)\n', N_parallel);
fprintf('    Actual total cells:      %d\n', N_cells_real);
fprintf('    Effective array area:    %.2f m²\n\n', A_sa_real);

% --- 3.6  Power generated by the sized array ---
P_sa_bol = A_sa_real * P0 * epsilon * Id;           % BoL power [W]
P_sa_eol = P_sa_bol * (1 - dpy)^mission_yrs;       % EoL power [W]

fprintf('  Power from sized solar array:\n');
fprintf('    BoL power (perp. to Sun): %.0f W\n', P_sa_bol);
fprintf('    EoL power:                %.0f W\n\n', P_sa_eol);

%% =========================================================================
%  SECTION 4: BATTERY SIZING
%  =========================================================================

fprintf('=== BATTERY SIZING ===\n');

% --- 4.1  Required energy capacity (Eq. 12) ---
%
%          Pe * Te / Xe
%  Ebatt = --------------
%           N * eta * DoD

N_packs = 2;    % Number of battery packs [-]
eta_bat  = 0.80; % Battery EoL efficiency [-]
DoD      = 0.30; % Depth of discharge [-]

E_batt = (P_eclipse * T_ecl_h / X_e) / (N_packs * eta_bat * DoD);  % [Wh]

fprintf('  Required battery energy capacity per pack (Eq. 12):\n');
fprintf('    Pe = %.0f W,  Te = %.2f h,  Xe = %.2f\n', P_eclipse, T_ecl_h, X_e);
fprintf('    N = %d packs,  eta = %.2f,  DoD = %.2f\n', N_packs, eta_bat, DoD);
fprintf('    Ebatt = (Pe*Te/Xe) / (N*eta*DoD) = %.0f Wh\n\n', E_batt);

% --- 4.2  Cell sizing (reference Li-Ion cell) ---
C_cell    = 3.0;    % Cell capacity [Ah]
V_cell_b  = 3.6;    % Cell voltage [V]
mu        = 0.80;   % Packing efficiency [-]

% Cells in series per string (to match bus voltage)
N_series_b = ceil(V_bus_min / V_cell_b);

% Capacity and energy of a single series string
C_string = mu * C_cell;                          % [Ah]
E_string = C_string * (V_cell_b * N_series_b);  % [Wh]

% Number of parallel strings (+ 1 redundant)
N_parallel_b = ceil(E_batt / E_string) + 1;

% Total cell count and actual energy/charge
N_cells_b    = N_parallel_b * N_series_b;
E_batt_real  = N_parallel_b * E_string;         % [Wh] per pack
C_batt       = N_parallel_b * mu * C_cell;      % [Ah] per pack

fprintf('  Cell configuration (Li-Ion: %.1f Ah, %.1f V; packing eff. %.2f):\n', ...
    C_cell, V_cell_b, mu);
fprintf('    Cells in series:             %d  (Vbus/Vcell = %g/%.1f)\n', ...
    N_series_b, V_bus_min, V_cell_b);
fprintf('    String capacity:             %.1f Ah\n', C_string);
fprintf('    Energy per string:           %.2f Wh\n', E_string);
fprintf('    Parallel strings:            %d  (incl. 1 redundant string)\n', N_parallel_b);
fprintf('    Total cells per pack:        %d\n', N_cells_b);
fprintf('    Energy capacity (per pack):  %.0f Wh\n', E_batt_real);
fprintf('    Electric charge (per pack):  %.2f Ah\n\n', C_batt);

% Total across all packs
fprintf('  Total battery (all %d packs combined):\n', N_packs);
fprintf('    Total energy:   %.0f Wh\n', N_packs * E_batt_real);
fprintf('    Total charge:   %.2f Ah\n\n', N_packs * C_batt);

%% =========================================================================
%  SECTION 5: SUMMARY TABLE
%  =========================================================================

fprintf('=================================================================\n');
fprintf('                        EPS SIZING SUMMARY\n');
fprintf('=================================================================\n');
fprintf('SOLAR ARRAYS\n');
fprintf('  Required power at EoL          %8.1f  W\n',   P_sa);
fprintf('  Spec. power at BoL             %8.1f  W/m²\n', P_bol_spec);
fprintf('  Spec. power at EoL             %8.1f  W/m²\n', P_eol_spec);
fprintf('  Required area                  %8.2f  m²\n',  A_sa);
fprintf('  Cells in series / parallel     %8d / %d\n',   N_series, N_parallel);
fprintf('  Total cell count               %8d\n',        N_cells_real);
fprintf('  Effective area                 %8.2f  m²\n',  A_sa_real);
fprintf('  BoL power (perp.)              %8.0f  W\n',   P_sa_bol);
fprintf('  EoL power                      %8.0f  W\n',   P_sa_eol);
fprintf('-----------------------------------------------------------------\n');
fprintf('BATTERIES  (per pack, N = %d packs)\n', N_packs);
fprintf('  Required energy per pack       %8.0f  Wh\n',  E_batt);
fprintf('  Cells in series / parallel     %8d / %d\n',   N_series_b, N_parallel_b);
fprintf('  Total cells per pack           %8d\n',        N_cells_b);
fprintf('  Actual energy per pack         %8.0f  Wh\n',  E_batt_real);
fprintf('  Electric charge per pack       %8.2f  Ah\n',  C_batt);
fprintf('  Total system energy            %8.0f  Wh\n',  N_packs * E_batt_real);
fprintf('  Total system charge            %8.2f  Ah\n',  N_packs * C_batt);
fprintf('=================================================================\n');
