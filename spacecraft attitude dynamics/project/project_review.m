clear; close all; clc;

%% satellite parameters
% inertia
ix = 12.1057;                     % [kg m^2]
iy = 15.0144;                     % [kg m^2]
iz = 18.0575;                     % [kg m^2]
I = diag([ix iy iz]);             % [kg m^2]

% DCM
A_BN0 = eye(3);

%% orbit parameters
% Keplerian parameters
a = 6900;                         % SMA [km]
e = 0.01;                         % eccentricity [-]
inc = deg2rad(10);                % inclination [rad]
TA0 = 0;                          % initial true anomaly [rad]

M_E = 5.97e24;                    % earth mass [kg]
G_E = 6.67e-20;                   % gravitational constant [km^3 kg^-1 s-^2]
mu_E = G_E*M_E;                   % earth gravitational parameter [km^3 s^-2]
R_E = 6378;                       % earth radius [km]
w_E = 7.27 * 10^-5;               % earth rotation rate [rad s^-1]
n = sqrt(mu_E / (a^3));           % mean rotational rate [rad s^-1]
T = 2*pi / n;                     % orbital period [s]

% angular velocities
w0 = [0; 0; n;]; % [rad s^-1]

%% sensors - gyroscope
% misalignment matrix
sens_gyro_theta_eps_x=1e-6;                 % [rad]
sens_gyro_theta_eps_y=1e-6;                 % [rad]
sens_gyro_theta_eps_z=1e-6;                 % [rad]

sens_gyro_theta_eps=[0,            -sens_gyro_theta_eps_z,  sens_gyro_theta_eps_y;
           sens_gyro_theta_eps_z,  0,             -sens_gyro_theta_eps_x;
           -sens_gyro_theta_eps_y, sens_gyro_theta_eps_x,    0];

% non-orthogonality
sens_gyro_eps_xy = 1e-6;  
sens_gyro_eps_xz = 1e-6; 
sens_gyro_eps_yz = 1e-6;  
sens_gyro_O = [1,      sens_gyro_eps_xy, sens_gyro_eps_xz;
     sens_gyro_eps_xy, 1,      sens_gyro_eps_yz;
     sens_gyro_eps_xz, sens_gyro_eps_yz, 1];

% dimensions
sens_gyro_m      = 0.052; % mass [kg]
sens_gyro_length = 44.8;  % length [mm] 
sens_gyro_height = 38.6;  % height [mm]
sens_gyro_depth  = 21.5;  % width [mm]

% datasheet values
sens_gyro_mis_err          = 1 / 1000;                     % misalignment error [mrad]
sens_gyro_run_run_bias     = 4 / 3600;                     % "sensor turn on" bias [deg s^-1]
sens_gyro_static_temp_bias = 9 / 3600;                     % "static temperature" bias [deg s^-1]
sens_gyro_SFE              = 500 * 1e-6;                   % scale factor [-]
sens_gyro_SFN              = 15 * 1e-6;                    % non linearity scale factor [-]
sens_gyro_update_rate      = 100;                          % [Hz]
sens_gyro_Ts               = 1 / sens_gyro_update_rate;    % sampling time [s]
sens_gyro_D_bias_inst      = 0.3 / 3600;                   % bias instability [deg s^-1]
sens_gyro_D_ARW            = 0.15 / sqrt(3600);            % angle random walk [deg s^-1/2]
sens_gyro_RRW              = 1e-3;                         % rate random walk
sens_gyro_D_RW             = sens_gyro_RRW / sqrt(3600);   % [deg s^-1]
sens_gyro_D_wn             = sens_gyro_D_ARW / sqrt(3600); % Amplitude of white noise 
sens_gyro_c_time           = 200;                          % [s]
sens_gyro_FS               = 10;                           % full scale [V]
sens_gyro_nbits            = 24;                           % [-] # bits
sens_gyro_LSB = sens_gyro_FS / (2 * exp(sens_gyro_nbits)); % Least significant bit

%% sensors - magnetometer
% misalignment matrix
sens_m_theta_eps_x=1e-6;                 % [rad]
sens_m_theta_eps_y=1e-6;                 % [rad]
sens_m_theta_eps_z=1e-6;                 % [rad]

sens_m_theta_eps=[0,            -sens_m_theta_eps_z,  sens_m_theta_eps_y;
           sens_m_theta_eps_z,  0,             -sens_m_theta_eps_x;
           -sens_m_theta_eps_y, sens_m_theta_eps_x,    0];

% non-orthogonality
sens_m_eps_xy = 1e-6;  
sens_m_eps_xz = 1e-6; 
sens_m_eps_yz = 1e-6;  
sens_m_O = [1,      sens_m_eps_xy, sens_m_eps_xz;
     sens_m_eps_xy, 1,      sens_m_eps_yz;
     sens_m_eps_xz, sens_m_eps_yz, 1];

% % dimensions
% sens_m_m=0.052;                  % [kg] mass of the gyroscope
% sens_m_length=44.8;              % [mm] 
% sens_m_height=38.6;              % [mm]
% sens_m_depth=21.5;               % [mm]
% 
% % datasheet values
% sens_m_mis_err=1/1000;           % misalignment error [mrad]
% sens_m_run_run_bias=4/3600;      % bias due to sensor turn on [deg s^-1]
% sens_m_static_temp_bias=9/3600;  % bias due to static temperature [deg s^-1]
% sens_m_SFE = 500*1e-6;           % scale factor [-]
% sens_m_SFN = 15*1e-6;            % non linearity scale factor [-]
% sens_m_update_rate = 100;        % [Hz]
% sens_m_Ts= 1/sens_m_update_rate;        % sampling time [s]
% sens_m_D_bias_inst=0.3/3600;     % bias instability [deg s^-1]
% sens_m_D_ARW=0.15/sqrt(3600);    % angle random walk [deg s^-1/2]
% sens_m_RRW=1e-3;                 % rate random walk GUESS
% sens_m_D_RW=sens_m_RRW/sqrt(3600);      % [deg s^-1]
% sens_m_D_wn = sens_m_D_ARW/sqrt(3600);  % Amplitude of white noise 
% sens_m_c_time=200;               % [s]
% sens_m_FS=10;                    % full scale [V]
% sens_m_nbits=24;                 % [-] # bits
% sens_m_LSB=sens_m_FS/(2*exp(sens_m_nbits));    % Least significant bit

%% sensors - earth horizon
% misalignment matrix
sens_eh_theta_eps_x = 1e-6; % [rad]
sens_eh_theta_eps_y = 1e-6; % [rad]
sens_eh_theta_eps_z = 1e-6; % [rad]

sens_eh_theta_eps=[0,            -sens_eh_theta_eps_z,  sens_eh_theta_eps_y;
           sens_eh_theta_eps_z,  0,             -sens_eh_theta_eps_x;
           -sens_eh_theta_eps_y, sens_eh_theta_eps_x,    0];

% non-orthogonality
sens_eh_eps_xy = 1e-6;  
sens_eh_eps_xz = 1e-6; 
sens_eh_eps_yz = 1e-6;  
sens_eh_O = [1,      sens_eh_eps_xy, sens_eh_eps_xz;
     sens_eh_eps_xy, 1,      sens_eh_eps_yz;
     sens_eh_eps_xz, sens_eh_eps_yz, 1];

% dimensions
sens_eh_m      = 18 * 10^-3; % mass [kg]
sens_eh_length = 35;         % length [mm] 
sens_eh_height = 24;         % height [mm]
sens_eh_depth  = 20;         % width [mm]

% datasheet values
sens_eh_accuracy         = 0.997;                     % accuracy (3-sigma)
% sens_eh_mis_err          = 1 / 1000;                   % misalignment error [mrad]
% sens_eh_run_run_bias     = 4 / 3600;                   % "sensor turn on" bias [deg s^-1]
% sens_eh_static_temp_bias = 9 / 3600;                   % "static temperature" bias [deg s^-1]
% sens_eh_SFE              = 500 * 1e-6;                 % scale factor [-]
% sens_eh_SFN              = 15 * 1e-6;                  % non linearity scale factor [-]
sens_eh_update_rate      = 2;                          % [Hz]
% sens_eh_Ts               = 1 / sens_eh_update_rate;    % sampling time [s]
% sens_eh_D_bias_inst      = 0.3 / 3600;                 % bias instability [deg s^-1]
% sens_eh_D_ARW            = 0.15 / sqrt(3600);          % angle random walk [deg s^-1/2]
% sens_eh_RRW              = 1e-3;                       % rate random walk
% sens_eh_D_RW             = sens_eh_RRW / sqrt(3600);   % [deg s^-1]
% sens_eh_D_wn             = sens_eh_D_ARW / sqrt(3600); % Amplitude of white noise 
% sens_eh_c_time           = 200;                        % [s]
% sens_eh_FS               = 10;                         % full scale [V]
% sens_eh_nbits            = 24;                         % [-] # bits
% sens_eh_LSB = sens_eh_FS / (2 * exp(sens_eh_nbits));   % Least significant bit

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

%% disturbances - gravity gradient
% see simulink model

%% simulation
% simulation options
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.01";
sim_options.StartTime = "0";
sim_options.StopTime = num2str(round(T, 0));

% simulation outputs
disp("running sim")
simout = sim("project_review_simulink.slx", sim_options);
disp("sim complete")

% time
t = simout.tout;

% gravity gradient GG
T_GG = simout.T_GG.Data;

% magnetism M
b_N = simout.b_N.Data; % inertial magnetic flux density
T_M = simout.T_M.Data; % magnetic torque

