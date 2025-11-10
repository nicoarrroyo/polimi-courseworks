clear; close all; clc;

%% initialisation
% inertia
ix = 25;
iy = 35;
iz = 50;
I = diag([ ix iy iz ]); % kg m^2

% DCM
A_BN0 = eye(3);

% kepler
a = 7093; % SMA [km]
e = 0.00277; % eccentricity [-]
TA0 = deg2rad(0); % initial true anomaly [rad]
inc = deg2rad(20); % inclination [rad]

% orbiting celestial body (earth)
M_E = 5.9722 * 10^24; % earth mass [kg]
G_E = 6.67 * 10^-11 * 10^-9; % gravitational constant [km^3 kg^-1 s^-2]
mu = M_E * G_E; % earth gravitational parameter [km^3 s^-2]
n = sqrt(mu/(a^3)); % average rotational rate [rad s^-1]
T = 2*pi / n; % period [s]

% angular velocities
w0 = [1e-6; 1e-6; n;];

% sun
n_sun = 2*pi / (60*60*24*365.25); % average rotation rate (sun) [rad s^-1]
r_sun = 1.496e8; % earth orbit radius (sun) [km]
e_sun = deg2rad(23.45); % eccentricity (sun) [-]

% SRP
F_e = 1358; % solar radiation intensity [W m^-2]
c = 2.998 * 10^8; % speed of light [m s^-1]
P = F_e / c;

% physical satellite properties
N_B = [
    [1 0 0] % body surface normal vectors
    [0 1 0]
    [-1 0 0]
    [0 -1 0] % unknown / assumed to follow 3
    [0 0 1] % unknown / assumed to follow 6
    [0 0 -1]
    [1 0 0] % solar panel surface normal vectors
    [-1 0 0]
    [1 0 0]
    [-1 0 0] % unknown / assumed to follow 3
    ];
rho_s = [0.5 0.5 0.5 0.5 0.5 0.5 ... % body surface solar values "s"
    0.1 0.1 0.1 0.1]; % solar panel surface solar values "s"
rho_d = [0.1 0.1 0.1 0.1 0.1 0.1 ... % body surface solar values "d"
    0.1 0.1 0.1 0.1]; % solar panel surface solar values "d"
area = [
    6; 6; 6; 6; 4; 4; ... % body surface areas [m^2]
    12; 12; 12; 12; % solar panel surface areas [m^2]
    ] * 10^-2;
r_Fi = [
    [10 0 0] % distances from body surface CoM to satellite CoM
    [0 10 0]
    [-10 0 0]
    [0 -10 0] % unknown / assumed to follow 3
    [0 0 15] % unknown / assumed to follow 6
    [0 0 -15]
    [0 45 0] % distances from body surface CoM to satellite CoM
    [0 45 0]
    [0 -45 0]
    [0 -45 0]
    ] * 10^-2;

% magnetism
j = [0.01; 0.05; 0.01;]; % magnetic dipole moment [amp m^2]
g_10 = -29404.8; % DGRF order 1 coefficients [nano tesla nT]
g_11 = -1450.9;
h_11 = 4652.5;
R_earth = 6378;
w_earth = 7.27 * 10^-5;

% gravity gradient estimate
I_M = max(max(I)); % maximum inertia moment
I_m = min(max(I)); % minimum inertia moment
max_GG = 3 * G_E * M_E / (2 * a^3) * abs(I_M - I_m);

% solar radiation pressure estimate
% note - this calculates a "worst-case" maximum everything
q = (1 + max(rho_s)) + (2 / 3) * max(rho_d) - 1;
max_arm = max(max(r_Fi));
max_SRP = P * max(area) * (1 + q) * (max_arm);

% magnetic torque estimate
H_0 = sqrt(g_10^2 + g_11^2 + h_11^2);
r_p = a * (1 - e);
B_max = (2 * (R_earth^3) * H_0) / (r_p^3) * 10^-9;
max_M = norm(j) * B_max;

% aerodynamic torque estimate
rho = 2.82e-14; % [kg/m^3] - ESTIMATE at 700km
C_D = 2.2; % [-] ASSUMPTION for LEO
A_s = 6 * 10^-2; % [m^2] - ASSUMPTION
lever_arm = 0.05; % 5 cm offset - ASSUMPTION
V = sqrt((mu * 10^9) / a * 10^3); % [m/s]
max_drag = 0.5 * rho * V^2 * A_s * C_D; % [N]
max_aero = max_drag * lever_arm;

% estimates outputs
fprintf("Estimated T_max (GG):   %.2e Nm\n", max_GG);
fprintf("Estimated T_max (SRP):  %.2e Nm\n", max_SRP);
fprintf("Estimated T_max (M):    %.2e Nm\n", max_M);
fprintf("Estimated T_max (aero): %.2e Nm\n", max_aero);

% simulation options
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.1";
sim_options.StartTime = "0";
sim_options.StopTime = "10";

%% outputs
disp("running sim")
simout = sim("project_review_simulink.slx", sim_options);
disp("sim complete")

% time
t = simout.tout;

% gravity gradient GG
T_GG = simout.T_GG.Data;

% solar radiation pressure SRP
S_N = simout.S_N.Data;
S_B = simout.S_B.Data;
S_B_hat = simout.S_B_hat.Data;
F_i = simout.F_i.Data; % SRP force
T_SRP = simout.T_i.Data; % SRP torque

% magnetism M
b_N = simout.b_N.Data; % inertial magnetic flux density
T_M = simout.T_M.Data; % magnetic torque

% total
T_tot = simout.T_tot.Data; % total torque
