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

T_orb   = 2*pi*sqrt(a^3/mu_E);% Orbit period            [sec]
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
fprintf(' Eclipse eff. Xe:      %.2f\n', X_ecl);
fprintf(' Sunlight eff. Xs:     %.2f\n\n', X_sun);

%% ========================================================================
%  SECTION 3: SOLAR ARRAY SIZING
%  ========================================================================

fprintf('======== SOLAR ARRAY SIZING ========\n');

% --- 3.1  Required power from solar arrays at EoL ---
sa.P = ((P_nam_ecl * T_ecl / X_ecl) + (P_nam_sun * T_sun / X_sun)) ...
       / T_sun;

% --- 3.2  Specific power at BoL ---
P0              = 1365;             % Solar const. at Earth         [W m^{-2}]
sa.epsilon      = 0.15;             % Silicon cell efficiency       [-]
sa.Id           = 0.80;             % Inherent degradation factor   [-]
sa.SAA_rad      = deg2rad(90-beta); % Sun Aspect Angle              [rad]

sa.P_bol_spec = P0 * sa.epsilon * sa.Id * cos(sa.SAA_rad); % Specific BoL power [W m^{-2}]

% --- 3.3  Specific power at EoL ---
sa.dpy          = 0.0275;           % Degradation per year          [-] (sl.34)
sa.mission_yrs  = 3;                % Mission duration              [years]

sa.P_eol_spec   = sa.P_bol_spec * (1 - sa.dpy)^sa.mission_yrs; % Specific EoL power  [W m^{-2}]

% --- 3.4  Required solar array area [m^2] ---
sa.A            = sa.P / sa.P_eol_spec;

% --- 3.5  Cell sizing (Silicon cells) ---
sa.cell_A       = 8*4e-4;           % Cell area                     [m^2] https://connectivity.esa.int/archives/projects/thin-film-solar-cell-demonstration-module
sa.cell_V       = 0.49;             % Cell max-power voltage        [V] https://ocw.tudelft.nl/wp-content/uploads/solar_energy_section_15_1-15_3.pdf, https://ntrs.nasa.gov/citations/19820061407, https://global.sharp/solar/en/space-qualified/pdf/datasheet_Si-CIC.pdf

temp_deg        = 2.2e-3;           % Voltage drop rate from temp.  [V/deg C]
temp_op         = 75;               % Operating temperature         [deg C]
temp_test       = 28;               % Temperature at cell test      [deg C]
temp_diff       = temp_op-temp_test;
temp_V_drop     = temp_deg*temp_diff;

sa.cell_V_real  = sa.cell_V - temp_V_drop;

% Minimum total cell and string count
sa.cell_Nmin    = ceil(sa.A / sa.cell_A);               % Min. total N. cells
sa.cell_N       = ceil(V_bus / sa.cell_V_real);         % N. series cells to match V
sa.strings_N    = ceil(sa.cell_Nmin / sa.cell_N) + 1;   % N. parallel strings (1 redundant)

% Actual total cell count and array area
sa.cell_N_real  = sa.strings_N * sa.cell_N;     % Real total N. cells
sa.A_real       = sa.cell_A * sa.cell_N_real;   % Real surface area [m^2]

% --- 3.6  Power generated by the sized array ---
sa.P_bol        = sa.A_real * P0 * sa.epsilon * sa.Id;   % BoL power [W]
sa.P_eol        = sa.P_bol * (1 - sa.dpy)^sa.mission_yrs;   % EoL power [W]

fprintf(' Required solar array power at EoL:\n');
fprintf('   P_sa    = ((Pe * Te / Xe) + (Ps * Ts / Xs)) / Ts\n');
fprintf('           = ((%.1f*%.1f/%.2f) + (%.1f*%.1f/%.2f)) / %.1f\n', ...
    P_nam_ecl, T_ecl, X_ecl, P_nam_sun, T_sun, X_sun, T_sun)
fprintf('           = %.1f W\n\n', sa.P);

fprintf(' Specific power at BoL:\n');
fprintf('   P_BoL   = P0 * eps * Id * cos(theta)\n');
fprintf('           = %.1f * %.2f * %.2f * cos(%.3f)\n', ...
    P0, sa.epsilon, sa.Id, sa.SAA_rad);
fprintf('           = %.1f W m^{-2}\n\n', sa.P_bol_spec)

fprintf(' Specific power at EoL:\n');
fprintf('   P_EoL   = P_BoL * (1 - dpy)^mission_yrs\n')
fprintf('           = %.1f * (1 - %.4f)^%.0f\n', sa.P_bol_spec, sa.dpy, sa.mission_yrs)
fprintf('           = %.1f W m^{-2}\n\n', sa.P_eol_spec)

fprintf(' Required solar array area:\n');
fprintf('   Asa     = P_sa / Peol\n');
fprintf('           = %.1f / %.1f\n', sa.P, sa.P_eol_spec);
fprintf('           = %.2f m^2\n\n', sa.A)

fprintf(' Cell configuration (silicon: A = %.2f cm², Vcell = %.2f V):\n', ...
    sa.cell_A*1e4, sa.cell_V_real);
fprintf('   Minimum cells required: %.0f\n', sa.cell_Nmin);
fprintf('   Cells in series:        %.0f  (Vbus/Vcell = %g/%.2f)\n', ...
    sa.cell_N, V_bus, sa.cell_V_real);
fprintf('   Strings in parallel:    %.0f  (incl. 1 redundant string)\n', sa.strings_N);
fprintf('   Actual total cells:     %.0f\n', sa.cell_N_real);
fprintf('   Effective array area:   %.2f m²\n\n', sa.A_real);

fprintf(' Power from sized solar array (perp. to Sun):\n');
fprintf('   BoL power:              %.0f W\n', sa.P_bol);
fprintf('   EoL power:              %.0f W\n\n', sa.P_eol);

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
ba.N_packs  = 1;    % N. battery packs      [-]
ba.eta      = 0.80; % Batt. EoL efficiency  [-] TODO
ba.DoD      = 0.15; % Depth of discharge    [-] cit. from franci

ba.E        = (P_nam_ecl * T_ecl / X_ecl) / (ba.N_packs * ba.eta * ba.DoD) / 3600; % [Wh]

% --- 4.2  Cell sizing (reference Li-Ion cell) ---
ba.cell_C       = 26.0; % Cell capacity         [Ah] 3 cells/package, 9 packages
ba.cell_V       = 3.3;  % Cell voltage          [V] up to 4.1
ba.cell_mu      = 0.80; % Packing efficiency    [-]

% Cell packages in series (to reach a reasonable bus voltage)
ba.package_N    = ceil(V_bus / ba.cell_V);

% Cells in parallel per package (to reach battery capacity)
ba.cell_C_real  = ba.cell_mu * ba.cell_C;
ba.package_E    = ba.cell_C_real * ba.package_N * ba.cell_V;
ba.cell_N       = ceil(ba.E / ba.package_E);

% Total cell count and actual energy and charge
ba.cell_Ntot    = ba.cell_N * ba.package_N;
ba.E_real       = ba.cell_N * ba.package_E;             % [Wh] per pack
ba.C_real       = ba.cell_N * ba.cell_mu * ba.cell_C;   % [Ah] per pack

fprintf(' Required battery energy capacity per pack:\n');
fprintf('   Ebatt   = (Pe*Te/Xe) / (N*eta*DoD) = %.0f Wh\n', ba.E);
fprintf('           = (%.1f*%.1f/%.1f) / (%.0f*%.2f*%.2f)\n', ...
    P_nam_ecl, T_ecl, X_ecl, ba.N_packs, ba.eta, ba.DoD);
fprintf('           = %.0f Wh\n\n', ba.E);

fprintf(' Cell configuration (Li-Ion: %.1f Ah, %.1f V; packing eff. %.2f):\n', ...
    ba.cell_C, ba.cell_V, ba.cell_mu);
fprintf('   Packages in series:     %.0f  (Vbus/Vcell = %g/%.1f)\n', ...
    ba.package_N, V_bus, ba.cell_V);
fprintf('   Cells per package:      %.0f\n',        ba.cell_N); %TODO 1 FOR REDUNDANCY?
fprintf('   Energy capacity:        %.0f Wh\n',     ba.E_real);
fprintf('   Electric charge:        %.2f Ah\n\n',   ba.C_real);

%% ========================================================================
%  SECTION 5: SUMMARY TABLE
%  ========================================================================
fprintf('=============================================================\n');
fprintf('                    EPS SIZING SUMMARY\n');
fprintf('=============================================================\n');
fprintf(' SOLAR ARRAYS\n');
fprintf('   Required power at EoL   %.1f  W\n',     sa.P);
fprintf('   Spec. power at BoL      %.1f  W/m²\n',  sa.P_bol_spec);
fprintf('   Spec. power at EoL      %.1f  W/m²\n',  sa.P_eol_spec);
fprintf('   Required area           %.2f  m²\n',    sa.A);
fprintf('   Cells in series         %.0f\n',        sa.cell_N);
fprintf('   Strings in parallel     %.0f\n',        sa.strings_N);
fprintf('   Total cell count        %.0f\n',        sa.cell_N_real);
fprintf('   Effective area          %.2f  m²\n',    sa.A_real);
fprintf('   BoL power (perp.)       %.0f  W\n',     sa.P_bol);
fprintf('   EoL power               %.0f  W\n',     sa.P_eol);
fprintf('-------------------------------------------------------------\n');
fprintf(' BATTERIES\n');
fprintf('   Required total energy   %.0f  Wh\n',    ba.E);
fprintf('   Cells in parallel       %.0f\n',        ba.cell_N);
fprintf('   Packages in series      %.0f\n',        ba.package_N);
fprintf('   Total cells             %.0f\n',        ba.cell_Ntot);
fprintf('   Actual total energy     %.f  Wh\n',     ba.E_real);
fprintf('   Total electric charge   %.2f  Ah\n',    ba.C_real);
fprintf('=============================================================\n');
