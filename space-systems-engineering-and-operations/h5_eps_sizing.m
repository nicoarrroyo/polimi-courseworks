%% EPS Sizing
% This script sizes the Electrical Power System (EPS) of the SMOS 
% spacecraft in a Dusk-Dawn Sun-Synchronous Orbit (SSO), covering:
%   1. Solar array sizing (area, cell count, series/parallel configuration)
%   2. Battery sizing (energy capacity, cell count, series/parallel config)

clc; clear; close all;

%% ========================================================================
%  SECTION 0: MISSION & ORBIT PARAMETERS
%  ========================================================================

R_E     = 6371;                 % Earth radius          [km]
tilt_E  = 23.44;                % Earth axial tilt      [deg]
mu_E    = 398600.4;             % Earth grav. param.    [km^2 s^{-3}]

h       = 765;                  % Altitude              [km]
a       = h + R_E;              % Semi-major axis       [km]
e       = 0.0001;               % Eccentricity          [-] (negligible)
i       = 98.44;                % Inclination           [deg]

T_orb   = 2*pi*sqrt(a^3/mu_E);% Orbit period          [sec]
T_orb_m = T_orb / 60;           % Orbital period        [min]
T_orb_h = T_orb_m / 60;         % Orbital period        [h]

beta    = 180 - tilt_E - i;

% Nominal orbit case is not considered for sizing

% Worst eclipse case (around winter)
T_ecl   = (T_orb / pi) * acos(sqrt(a^2 - R_E^2) / (a * cosd(beta)));

T_sun   = T_orb - T_ecl;

fprintf('======== ORBIT PARAMETERS ========\n');
fprintf(' Semi-Major Axis:      %.2f km\n', a);
fprintf(' Orbital period:       %.2f minutes\n', T_orb/60);
fprintf(' Winter sunlit time:   %.2f minutes\n', T_sun/60);
fprintf(' Winter eclipse time:  %.2f minutes\n\n', T_ecl/60);

%% ========================================================================
%  SECTION 1: POWER BUDGET (from Table 2)
%  ========================================================================
% Safe Hold, Star Acquisition, Normal Autonomous, Orbit Correction 2/4

% Safe-Hold Mode [W]
P_shm       = 110*1.2;

% Star Acquisition Mode [W]
P_sam       = 248*1.2;

% Normal Autonomous Mode nominally [W]
P_nam_sun   = 665*1.2;

% Normal Autonomous Mode during eclipse [W]
P_nam_ecl   = 765*1.2; % TODO heaters at 160 W instead of 80 W

% Orbit COrrection Mode 2 [W]
P_ocm2      = 646.4*1.2;

% Orbit Correction Mode 4 [W]
P_ocm4      = 686.3*1.2;

fprintf('======== POWER BUDGET (w/ 20%% margin) ========\n');
fprintf(' SHM:                  %.1f W\n', P_shm);
fprintf(' SAM:                  %.1f W\n', P_sam);
fprintf(' NAM (Sun.):           %.1f W\n', P_nam_sun);
fprintf(' NAM (Ecl.):           %.1f W\n', P_nam_ecl);
fprintf(' OCM2:                 %.1f W\n', P_ocm2);
fprintf(' OCM4:                 %.1f W\n\n', P_ocm4);

%% =========================================================================
%  SECTION 2: EPS CONFIGURATION & EFFICIENCIES
%  =========================================================================

% Unregulated bus voltages min and max [V]
V_bus_min   = 23; % TODO
V_bus_max   = 37;
V_bus       = 28; % TODO reasonable guess??

% Line efficiencies DET
X_ecl       = 0.65;
X_sun       = 0.85;

fprintf('======== EPS CONFIGURATION ========\n');
fprintf(' Eclipse efficiency Xe:    %.2f\n', X_ecl);
fprintf(' Sunlight efficiency Xs:   %.2f\n\n', X_sun);

%% ========================================================================
%  SECTION 3: SOLAR ARRAY SIZING
%  ========================================================================

fprintf('======== SOLAR ARRAY SIZING ========\n');

% --- 3.1  Required power from solar arrays at EoL ---
P_sa = ((P_nam_ecl * T_ecl / X_ecl) + (P_nam_sun * T_sun / X_sun)) ...
       / T_sun;

% --- 3.2  Specific power at BoL ---
P0          = 1365;             % Solar const. at Earth         [W m^{-2}]
% epsilon 0.25 https://www.nlr.gov/pv/cell-efficiency
epsilon     = 0.15;             % Silicon cell efficiency       [-]
Id          = 0.80;             % Inherent degradation factor   [-]
SAA_rad     = deg2rad(90-beta); % Sun Aspect Angle              [rad]

P_bol_spec = P0 * epsilon * Id * cos(SAA_rad); % Specific BoL power [W m^{-2}]

% --- 3.3  Specific power at EoL ---
dpy             = 0.0275;       % Degradation per year          [-] (sl.34)
mission_yrs     = 3;            % Mission duration              [years]

P_eol_spec      = P_bol_spec * (1 - dpy)^mission_yrs; % Specific EoL power  [W m^{-2}]

% --- 3.4  Required solar array area ---
A_sa            = P_sa / P_eol_spec;    % Required solar array area         [m^2]

% --- 3.5  Cell sizing (Silicon cells) ---
A_cell_sa       = 8*4e-4;       % Cell area                 [m^2] https://connectivity.esa.int/archives/projects/thin-film-solar-cell-demonstration-module
V_cell_sa       = 0.55;         % Cell max-power voltage    [V] https://ocw.tudelft.nl/wp-content/uploads/solar_energy_section_15_1-15_3.pdf, https://ntrs.nasa.gov/citations/19820061407, https://global.sharp/solar/en/space-qualified/pdf/datasheet_Si-CIC.pdf

% Minimum total cell and string count
N_cells_min_sa  = ceil(A_sa / A_cell_sa);   % Min. total N. cells
N_series_sa     = ceil(V_bus / V_cell_sa);   % N. series cells to match V
N_parallel_sa   = ceil(N_cells_min_sa / N_series_sa) + 1; % N. parallel strings (1 redundant)

% Actual total cell count and array area
N_cells_real_sa = N_parallel_sa * N_series_sa;  % Real total N. cells
A_sa_real       = A_cell_sa * N_cells_real_sa;  % [m^2]

% --- 3.6  Power generated by the sized array ---
P_sa_bol = A_sa_real * P0 * epsilon * Id;       % BoL power [W]
P_sa_eol = P_sa_bol * (1 - dpy)^mission_yrs;    % EoL power [W]

fprintf(' Required solar array power at EoL:\n');
fprintf('   P_sa    = ((Pe * Te / Xe) + (Ps * Ts / Xs)) / Ts\n');
fprintf('           = ((%.1f*%.1f/%.2f) + (%.1f*%.1f/%.2f)) / %.1f\n', ...
    P_nam_ecl, T_ecl, X_ecl, P_nam_sun, T_sun, X_sun, T_sun)
fprintf('           = %.1f W\n\n', P_sa);

fprintf(' Specific power at BoL:\n');
fprintf('   P_BoL   = P0 * eps * Id * cos(theta)\n');
fprintf('           = %.1f * %.1f * %.1f * cos(%.3f)\n', ...
    P0, epsilon, Id, SAA_rad);
fprintf('           = %.1f W m^{-2}\n\n', P_bol_spec)

fprintf(' Specific power at EoL:\n');
fprintf('   P_EoL   = P_BoL * (1 - dpy)^mission_yrs\n')
fprintf('           = %.1f * (1 - %.2f)^%.0f\n', P_bol_spec, dpy, mission_yrs)
fprintf('           = %.1f W m^{-2}\n\n', P_eol_spec)

fprintf(' Required solar array area:\n');
fprintf('   Asa     = P_sa / Peol\n');
fprintf('           = %.1f / %.1f\n', P_sa, P_eol_spec);
fprintf('           = %.2f m^2\n\n', A_sa)

fprintf(' Cell configuration (silicon: A = %.2f cm², Vcell = %.2f V):\n', ...
    A_cell_sa*1e4, V_cell_sa);
fprintf('   Minimum cells required:  %d\n', N_cells_min_sa);
fprintf('   Cells in series:         %d  (Vbus/Vcell = %g/%.2f)\n', ...
    N_series_sa, V_bus, V_cell_sa);
fprintf('   Parallel strings:        %d  (incl. 1 redundant string)\n', N_parallel_sa);
fprintf('   Actual total cells:      %d\n', N_cells_real_sa);
fprintf('   Effective array area:    %.2f m²\n\n', A_sa_real);

fprintf(' Power from sized solar array:\n');
fprintf('   BoL power (perp. to Sun): %.0f W\n', P_sa_bol);
fprintf('   EoL power:                %.0f W\n\n', P_sa_eol);

%% ========================================================================
%  SECTION 4: BATTERY SIZING
%  ========================================================================
% SMOS uses a Parallel-Series battery architecture. 3 cells in parallel
% with each other, which makes up a cell package. Then there are 9 cell
% packages, all in series with each other. These 9 cell packages (a total
% of 27 cells) make up the battery. eddajje
V_min = 24;
V_max = 37;
V_cell_min = 3.3;
V_cell_max = 4.1;
% nickel cadmium?
% https://www.spiedigitallibrary.org/conference-proceedings-of-spie/10570/2326414/A-new-European-small-platform--Proteus-and-prospected-optical/10.1117/12.2326414.full

fprintf('======== BATTERY SIZING ========\n');

% --- 4.1  Required energy capacity ---
% shall we use a lithium thionyl chloride cell? page 13 of slides
N_packs         = 1;    % N. battery packs      [-]
eta_bat         = 0.80; % Batt. EoL efficiency  [-] TODO
DoD             = 0.30; % Depth of discharge    [-] TODO

E_batt = (P_nam_ecl * T_ecl / X_ecl) / (N_packs * eta_bat * DoD) / 3600; % [Wh]

fprintf(' Required battery energy capacity per pack:\n');
fprintf('   Ebatt   = (Pe*Te/Xe) / (N*eta*DoD) = %.0f Wh\n', E_batt);
fprintf('           = (%.1f*%.1f/%.1f) / (%.0f*%.2f*%.2f)\n', ...
    P_nam_ecl, T_ecl, X_ecl, N_packs, eta_bat, DoD);
fprintf('           = %.0f Wh\n\n', E_batt);

% --- 4.2  Cell sizing (reference Li-Ion cell) ---
C_cell_batt     = 26.0; % Cell capacity         [Ah] 3 cells/package, 9 packages
V_cell_batt     = 3.3;  % Cell voltage          [V] up to 4.1
pack_batt       = 0.80; % Packing efficiency    [-]

% Cell packages in series (to reach a reasonable bus voltage)
N_series_batt   = ceil(V_bus / V_cell_batt);

% Cells in parallel per package (to reach battery capacity)
C_package       = pack_batt * C_cell_batt;
E_package       = C_package * N_series_batt * V_cell_batt;
N_parallel_batt = ceil(E_batt / E_package);

% Total cell count and actual energy/charge
N_cells_batt    = N_parallel_batt * N_series_batt;
E_batt_real     = N_parallel_batt * E_package;                  % [Wh] per pack
C_batt          = N_parallel_batt * pack_batt * C_cell_batt;    % [Ah] per pack

fprintf(' Cell configuration (Li-Ion: %.1f Ah, %.1f V; packing eff. %.2f):\n', ...
    C_cell_batt, V_cell_batt, pack_batt);
fprintf('   Packages in series:          %d  (Vbus/Vcell = %g/%.1f)\n', ...
    N_series_batt, V_bus, V_cell_batt);
fprintf('   Package capacity:            %.1f Ah\n', C_package);
fprintf('   Energy per package:          %.2f Wh\n', E_package);
fprintf('   Total cells per pack:        %d\n', N_parallel_batt); %TODO 1 FOR REDUNDANCY?
fprintf('   Energy capacity (per pack):  %.0f Wh\n', E_batt_real);
fprintf('   Electric charge (per pack):  %.2f Ah\n\n', C_batt);

% Total across all packs
fprintf(' Total battery (all %d packs combined):\n', N_packs);
fprintf('   Total energy:   %.0f Wh\n', N_packs * E_batt_real);
fprintf('   Total charge:   %.2f Ah\n\n', N_packs * C_batt);

%% ========================================================================
%  SECTION 5: SUMMARY TABLE
%  ========================================================================
fprintf('=============================================================\n');
fprintf('                    EPS SIZING SUMMARY\n');
fprintf('=============================================================\n');
fprintf(' SOLAR ARRAYS\n');
fprintf('   Required power at EoL       %8.1f  W\n',   P_sa);
fprintf('   Spec. power at BoL          %8.1f  W/m²\n', P_bol_spec);
fprintf('   Spec. power at EoL          %8.1f  W/m²\n', P_eol_spec);
fprintf('   Required area               %8.2f  m²\n',  A_sa);
fprintf('   Cells in series / parallel  %8d / %d\n',   N_series_sa, N_parallel_sa);
fprintf('   Total cell count            %8d\n',        N_cells_real_sa);
fprintf('   Effective area              %8.2f  m²\n',  A_sa_real);
fprintf('   BoL power (perp.)           %8.0f  W\n',   P_sa_bol);
fprintf('   EoL power                   %8.0f  W\n',   P_sa_eol);
fprintf('-------------------------------------------------------------\n');
fprintf(' BATTERIES  (per pack, N = %d packs)\n', N_packs);
fprintf('   Required energy per pack    %8.0f  Wh\n',  E_batt);
fprintf('   Cells in series / parallel  %8d / %d\n',   N_series_batt, N_parallel_batt);
fprintf('   Total cells per pack        %8d\n',        N_cells_batt);
fprintf('   Actual energy per pack      %8.0f  Wh\n',  E_batt_real);
fprintf('   Electric charge per pack    %8.2f  Ah\n',  C_batt);
fprintf('   Total system energy         %8.0f  Wh\n',  N_packs * E_batt_real);
fprintf('   Total system charge         %8.2f  Ah\n',  N_packs * C_batt);
fprintf('=============================================================\n');
