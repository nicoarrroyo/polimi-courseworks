clear; close all; clc;

%% satellite parameters
% inertia
ix = 25;                          % [kg m^2]
iy = 35;                          % [kg m^2]
iz = 50;                          % [kg m^2]
I = diag([ix iy iz]);             % [kg m^2]

% DCM
A_BN0 = eye(3);

%% orbit parameters
% Keplerian parameters
a = 7000;                         % SMA [km]
e = 0.001;                        % eccentricity [-]
inc = deg2rad(20);                % inclination [rad]
TA0 = 0;                          % initial true anomaly [rad]

M_E = 5.97e24;                    % earth mass [kg]
G_E = 6.67e-20;                   % gravitational constant [km^3 kg^-1 s-^2]
mu_E = G_E*M_E;                   % earth gravitational parameter [km^3 s^-2]
R_E = 6378;                       % earth radius [km]
w_E = 7.27 * 10^-5;               % earth rotation rate [rad s^-1]
n = sqrt(mu_E/(a^3));             % mean rotational rate [rad s^-1]
T = 2*pi / n;                     % orbital period [s]

% angular velocities
w0 = [10^-6; 10^-6; n;]; % rad s^-1

%% sensors - gyroscope
% Initial misalignment angles
theta_eps_x=1e-6;                 % [rad]
theta_eps_y=1e-6;                 % [rad]
theta_eps_z=1e-6;                 % [rad]

theta_eps=[0,            -theta_eps_z,  theta_eps_y;
           theta_eps_z,  0,             -theta_eps_x;
           -theta_eps_y, theta_eps_x,    0];

% Non-Orthogonality matrix
eps_xy = 1e-6;  
eps_xz = 1e-6; 
eps_yz = 1e-6;  
O = [1,      eps_xy, eps_xz;
     eps_xy, 1,      eps_yz;
     eps_xz, eps_yz, 1];

% Data (from data sheet)
m=0.052;                  % [kg] mass of the gyroscope
length=44.8;              % [mm] 
height=38.6;              % [mm]
deep=21.5;                % [mm]

mis_err=1/1000;           % misalignment error [mrad]
run_run_bias=4/3600;      % bias due to sensor turn on [deg s^-1]
static_temp_bias=9/3600;  % bias due to static temperature [deg s^-1]
SFE = 500*1e-6;           % scale factor [-]
SFN = 15*1e-6;            % non linearity scale factor [-]
update_rate = 100;        % [Hz]
Ts= 1/update_rate;        % sampling time [s]
D_bias_inst=0.3/3600;     % bias instability [deg s^-1]
D_ARW=0.15/sqrt(3600);    % angle random walk [deg s^-1/2]
RRW=1e-3;                 % rate random walk GUESS
D_RW=RRW/sqrt(3600);      % [deg s^-1]
D_wn = D_ARW/sqrt(3600);  % Amplitude of white noise 
c_time=200;               % [s]
FS=10;                    % full scale [V]
nbits=24;                 % [-] # bits
LSB=FS/(2*exp(nbits));    % Least significant bit

%% disturbances - SRP
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

%% disturbances - magnetism
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

%% disturbances - magnetism order 5
% we chose order 5 because order 6 had less than a 1% difference in total 
% magnitude of the magnetic field vector (see mag_order_pathfinder.m)
lat = 45.4685; % current dummy latitude (over milan) [deg]
long = 9.1824; % current dummy longitude (over milan) [deg]
B_LVLH = get_mag_order_5(R_E, a, lat, long);

%% disturbance estimates
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
B_max = (2 * (R_E^3) * H_0) / (r_p^3) * 10^-9;
max_M = norm(j) * B_max;

% aerodynamic torque estimate
rho = 2.82e-14; % [kg/m^3] - ESTIMATE at 700km
C_D = 2.2; % [-] ASSUMPTION for LEO
A_s = 6 * 10^-2; % [m^2] - ASSUMPTION
lever_arm = 0.05; % 5 cm offset - ASSUMPTION
V = sqrt((mu_E * 10^9) / a * 10^3); % [m/s]
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
sim_options.FixedStep = "0.01";
sim_options.StartTime = "0";
sim_options.StopTime = num2str(round(T/5, 0));

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

