clear; close all; clc;

% inertia
ix = 19.5;
iy = 19;
iz = 12.6;
I = diag([ ix iy iz ]); % kg m^2

% DCM
A_BN0 = eye(3);

%% orbit parameters
% keplerian parameters
a = 630 + 6378; % SMA [km]
e = 0; % eccentricity [-]
TA0 = deg2rad(0); % initial true anomaly [rad]
inc = deg2rad(97.9); % inclination [rad]

% orbiting celestial body (earth)
M_E = 5.9722 * 10^24; % earth mass [kg]
G_E = 6.67 * 10^-11 * 10^-9; % gravitational constant [km^3 kg^-1 s^-2]
mu = M_E * G_E; % earth gravitational parameter [km^3 s^-2]
n = sqrt(mu/(a^3)); % average rotational rate [rad s^-1]
T = 2*pi / n; % period [s]

% angular velocities
w0 = [10^-6; 10^-6; n;]; % rad s^-1

%% disturbances
% sun
n_sun = 2*pi / (60*60*24*365.25); % average rotation rate (sun) [rad s^-1]
r_sun = 1.496e8; % earth orbit radius (sun) [km]
e_sun = deg2rad(23.45); % eccentricity (sun) [-]

% SRP
F_e = 1358; % solar radiation intensity [W m^-2]
c = 2.998 * 10^8; % speed of light [m s^-1]
P = F_e / c;

% physical satellite properties
N_B = [ % body surface normal vectors
    [1 0 0]
    [0 1 0]
    [-1 0 0]
    [0 -1 0]
    [0 0 1]
    [0 0 -1]
    ];
rho_s = [0.1 0.5 0.1 0.5 0.5 0.1]; % body surface solar values "s"
rho_d = [0.1 0.1 0.1 0.1 0.1 0.1]; % body surface solar values "d"
area = [ % body surface areas [m^2]
    1.08; 0.96; -1.08; -0.96; 0.72; -0.72;
    ] * 10^-2; % calculated from dimensions 0.8 m x 0.9 m x 1.2 m
r_Fi = [ % distances from body surface CoM to satellite CoM
    [120/2 0 0]
    [0 90/2 0]
    [-120/2 0 0]
    [0 -90 0]
    [0 0 80/2]
    [0 0 -80/2]
    ] * 10^-2; % assumes CoM is in the geometric centre of the satellite

% magnetism
j = [0.01; 0.05; 0.01;]; % magnetic dipole moment [amp m^2]
g_10 = -29404.8; % DGRF order 1 coefficients [nano tesla nT]
g_11 = -1450.9;
h_11 = 4652.5;
g_20 = -2499.6; % DGRF order 2 coefficients
g_21 = 2982.0;
g_22 = 1677.0;
h_21 = -2991.6;
h_22 = -734.6;
g_30 = 1363.2; % DGRF order 3 coefficients
g_31 = -2381.2;
g_32 = 1236.2;
g_33 = 525.7;
h_31 = -82.1;
h_32 = 241.9;
h_33 = -543.4;
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

%% sensors
% gyroscope (STIM210 Multi-Axis Gyro Module)
gyro_arw = 0.15; % angular random walk [deg s^-1 hr^(-1/2)]
gyro_arw = deg2rad(gyro_arw) / sqrt(3600); % [rad s^-1]

gyro_bias = 0.4; % static bias (bias instability) [deg hr^-1]
gyro_bias = deg2rad(gyro_bias) / 3600; % [rad s^-1]

gyro_misalignment = 1 * 10^-3; % misalignment (pg5) [rad]

% simulation options
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.1";
sim_options.StartTime = "0";
sim_options.StopTime = num2str(round(T, 0));

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

