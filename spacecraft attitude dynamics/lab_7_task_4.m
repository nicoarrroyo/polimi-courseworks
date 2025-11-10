clear; close all; clc;

%% initial conditions
% inertia
ix = 100.9;
iy = 25.1;
iz = 91.6;
I = diag([ ix iy iz ] * 10^-2); % kg m^2

% angular velocities
w0 = [0.45; 0.52; 0.55;];

% DCM
A_BN0 = eye(3);

% keplerian parameters
a = 7093; % SMA [km]
e = 0.00277; % eccentricity [-]
TA0 = deg2rad(0); % initial true anomaly [rad]
inc = deg2rad(98.27); % inclination [rad]

mu = astroConstants(13); % earth gravitational parameter [km^3 kg^-1 s^-2]
n = sqrt(mu/(a^3)); % average rotational rate [rad s^-1]
T = 2*pi / n; % period [s]

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
g_10 = -29404.8; % DGRF order 1 coefficients [nano tesla]
g_11 = -1450.9;
h_11 = 4652.5;
R_earth = 6378;
w_earth = 7.27 * 10^-5;

% simulation options
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.1";
sim_options.StartTime = "0";
sim_options.StopTime = "10";

%% outputs
disp("running sim")
simout = sim("lab_7_task_4_simulink", sim_options);
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
