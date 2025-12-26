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
Omega = 80.5;                     % Right Ascension of Ascending node [rad]
omega = 0;                        % Argument of Perigee [rad]
TA0 = 0;                          % initial true anomaly [rad]

M_E = 5.97e24;                    % earth mass [kg]
G_E = 6.67e-20;                   % gravitational constant [km^3 kg^-1 s-^2]
mu_E = G_E*M_E;                   % earth gravitational parameter [km^3 s^-2]
R_E = 6378;                       % earth radius [km]
w_E = 7.27 * 10^-5;               % earth rotation rate [rad s^-1]
n = sqrt(mu_E / (a^3));           % mean rotational rate [rad s^-1]
T = 2*pi / n;                     % orbital period [s]

% angular velocities
w0 = [0; 0; 0;]; % [rad s^-1]

%% Sensors - Gyroscope
% misalignment matrix
sens_gyro.theta_eps_x=1e-6;                 % [rad]
sens_gyro.theta_eps_y=1e-6;                 % [rad]
sens_gyro.theta_eps_z=1e-6;                 % [rad]

sens_gyro.theta_eps=[0,            -sens_gyro.theta_eps_z,  sens_gyro.theta_eps_y;
           sens_gyro.theta_eps_z,  0,             -sens_gyro.theta_eps_x;
           -sens_gyro.theta_eps_y, sens_gyro.theta_eps_x,    0];

% non-orthogonality
sens_gyro.eps_xy = 1e-6;  
sens_gyro.eps_xz = 1e-6; 
sens_gyro.eps_yz = 1e-6;  
sens_gyro.O = [1,      sens_gyro.eps_xy, sens_gyro.eps_xz;
     sens_gyro.eps_xy, 1,      sens_gyro.eps_yz;
     sens_gyro.eps_xz, sens_gyro.eps_yz, 1];

% DCM GYRO BODY (fixed matrix, doesnt change in time, in our case it is the
% identity matrix)
sens_gyro.dcm_body = eye(3);

% dimensions
sens_gyro.m      = 0.052; % mass [kg]
sens_gyro.length = 44.8;  % length [mm] 
sens_gyro.height = 38.6;  % height [mm]
sens_gyro.depth  = 21.5;  % width [mm]

% datasheet values
sens_gyro.mis_err          = 1 / 1000;                     % misalignment error [mrad]
sens_gyro.run_run_bias     = 4 / 3600;                     % "sensor turn on" bias [deg s^-1]
sens_gyro.static_temp_bias = 9 / 3600;                     % "static temperature" bias [deg s^-1]
sens_gyro.SFE              = 500 * 1e-6;                   % scale factor [-]
sens_gyro.SFN              = 15 * 1e-6;                    % non linearity scale factor [-]
sens_gyro.update_rate      = 100;                          % [Hz]
sens_gyro.Ts               = 1 / sens_gyro.update_rate;    % sampling time [s]
sens_gyro.D_bias_inst      = 0.3 / 3600;                   % bias instability [deg s^-1]
sens_gyro.D_ARW            = 0.15 / sqrt(3600);            % angle random walk [deg s^-1/2]
sens_gyro.RRW              = 1e-3;                         % rate random walk
sens_gyro.D_RW             = sens_gyro.RRW / sqrt(3600);   % [deg s^-1]
sens_gyro.D_wn             = sens_gyro.D_ARW / sqrt(3600); % Amplitude of white noise 
sens_gyro.c_time           = 200;                          % [s]
sens_gyro.FS               = 10;                           % full scale [V]
sens_gyro.nbits            = 24;                           % [-] # bits
sens_gyro.LSB = sens_gyro.FS / (2 ^ sens_gyro.nbits); % Least significant bit

%% Sensors - Magnetometer
% Based on AAC SpaceQuest MAG-3 3-axis fluxgate magnetometer
% Space-qualified, TRL-9, flight-proven on numerous missions

% --- Physical Dimensions ---
sens_mag.length = 82.6e-3;  % length [m]
sens_mag.width  = 35.1e-3;  % width [m]
sens_mag.height = 32.3e-3;  % height [m]
sens_mag.m      = 100e-3;   % mass [kg]

% --- Hard Iron and Soft Iron Effects ---
% Hard iron: constant magnetic field from spacecraft (ferromagnetic materials)
sens_mag.hard_iron = [50e-9; 30e-9; -40e-9];  % [T] typical spacecraft disturbance

% Soft iron: induced magnetization (scaling and rotation)
% Represented as a 3x3 matrix deviation from identity
sens_mag.soft_iron = eye(3) + 0.01 * randn(3);  % 1% soft iron effect

% --- Misalignment Matrix ---
% Installation misalignment (sensor to body frame)
sens_mag.theta_eps_x = deg2rad(0.5);  % [rad] typical mounting tolerance
sens_mag.theta_eps_y = deg2rad(0.5);  % [rad]
sens_mag.theta_eps_z = deg2rad(0.5);  % [rad]

sens_mag.theta_eps = [0,                     -sens_mag.theta_eps_z,  sens_mag.theta_eps_y;
                      sens_mag.theta_eps_z,   0,                    -sens_mag.theta_eps_x;
                     -sens_mag.theta_eps_y,   sens_mag.theta_eps_x,  0];

% Total misalignment transformation
sens_mag.C_misalign = eye(3) + sens_mag.theta_eps;

% --- Non-Orthogonality Matrix (between sensor axes) ---
% Based on ±1° orthogonality spec
sens_mag.eps_xy = deg2rad(0.5);  % [rad]
sens_mag.eps_xz = deg2rad(0.5);  % [rad]
sens_mag.eps_yz = deg2rad(0.5);  % [rad]

sens_mag.O = [1,                 sens_mag.eps_xy, sens_mag.eps_xz;
              sens_mag.eps_xy,   1,               sens_mag.eps_yz;
              sens_mag.eps_xz,   sens_mag.eps_yz, 1];

% --- Performance Specifications ---
% Accuracy (as percentage of full scale)
sens_mag.accuracy = 0.0075;  % ±0.75% of FS (0.5% typical)

% Field measurement range
sens_mag.range = 65e-6;  % ±65 μT (standard range) [T]

% Noise density
sens_mag.noise_density = 12e-12;  % 12 pT/√Hz @1Hz [T/√Hz]

% --- Temperature Effects ---
% Scale factor temperature shift
sens_mag.scale_factor_temp = 0.00007;  % 0.007% FS/°C

% Zero shift with temperature
sens_mag.zero_temp_shift = 0.6e-9;  % ±0.6 nT/°C [T/°C]

% --- Axial Alignment ---
% Orthogonality between axes
sens_mag.orthogonality_error = deg2rad(1);  % Better than ±1 deg [rad]

% --- Perming Effects (magnetic memory) ---
% Susceptibility to perming (residual magnetization)
sens_mag.perming_shift = 8e-9; % ±8 nT shift with ±5 Gauss applied [T]
sens_mag.perming_field = 5e-4; % ±5 Gauss = 500 μT [T]

% --- Frequency Response ---
sens_mag.bandwidth = 500; % 3dB @ >500 Hz (up to 4 kHz wideband)

% --- Sampling Rate ---
sens_mag.update_rate = 10; % [Hz] typical for attitude determination
sens_mag.Ts = 1 / sens_mag.update_rate; % sampling time [s]

% --- Digital Conversion (ADC) ---
sens_mag.nbits = 16; % 16-bit ADC typical for spacecraft magnetometers
sens_mag.LSB = (2 * sens_mag.range) / (2^sens_mag.nbits); % [T]

% --- Derived Specifications ---
% RMS noise at update rate
sens_mag.noise_rms = sens_mag.noise_density * sqrt(sens_mag.update_rate); % [T]

% 1-sigma measurement uncertainty (from accuracy spec)
sens_mag.sigma_1axis = sens_mag.accuracy * sens_mag.range; % [T]

% For Simulink Band-Limited White Noise block
sens_mag.noise_power = sens_mag.sigma_1axis^2; % [rad^2]

% Scale factor (converts measured field to true field)
sens_mag.scale_factor_error = 0.001; % nominal, will add error term

% --- Bias Instability (Random Walk) ---
% Not typically specified for magnetometers, but exists
sens_mag.bias_instability = 0.1e-9; % [T/s] estimated (0.1 nT/s)

% --- Output Flags ---
sens_mag.valid_threshold = 1.5 * sens_mag.range; % saturation threshold

%% Sensors - Earth Horizon
% Based on MAI-SES IR Earth Sensor specifications

% --- Physical Dimensions ---
sens_eh.m      = 33e-3;     % mass [kg]
sens_eh.length = 43.3e-3;   % length [m] - CORRECTED to meters
sens_eh.width  = 31.8e-3;   % width [m]
sens_eh.height = 20.7e-3;   % height [m]

% --- Optical Characteristics ---
sens_eh.fov_coarse       = deg2rad(60);      % coarse field of view [rad]
sens_eh.fov_fine         = deg2rad(7);       % fine field of view [rad]

% --- Accuracy Specifications (using fine FOV resolution) ---
sens_eh.resolution_fine   = 0.25;   % fine FOV resolution [deg]
sens_eh.resolution_coarse = 1.0;    % coarse FOV resolution [deg]

% Convert to radians for simulation
sens_eh.sigma_1axis = deg2rad(sens_eh.resolution_fine); % 1-sigma noise per axis [rad]

% For 2-axis measurement (pitch and roll from nadir)
sens_eh.sigma_pitch = sens_eh.sigma_1axis;  % [rad]
sens_eh.sigma_roll  = sens_eh.sigma_1axis;  % [rad]

% --- Bias and Systematic Errors ---
% Static bias (offset from true nadir direction)
sens_eh.bias_pitch = deg2rad(0.1);   % pitch bias [rad] (~0.1 deg typical)
sens_eh.bias_roll  = deg2rad(0.1);   % roll bias [rad]

% Bias instability (slow drift over time)
sens_eh.bias_instability = deg2rad(0.05/3600); % [rad/s] (0.05 deg/hr typical)

% --- Misalignment Errors ---
% Installation misalignment (sensor to body frame)
sens_eh.theta_eps_x = deg2rad(0.01);  % [rad] (~0.01 deg typical mounting tolerance)
sens_eh.theta_eps_y = deg2rad(0.01);  % [rad]
sens_eh.theta_eps_z = deg2rad(0.01);  % [rad]

% Misalignment matrix (small angle approximation)
% this is applicable because we won't use earth-horizon for detumbling
sens_eh.theta_eps = [0,                    -sens_eh.theta_eps_z,  sens_eh.theta_eps_y;
                     sens_eh.theta_eps_z,   0,                   -sens_eh.theta_eps_x;
                    -sens_eh.theta_eps_y,   sens_eh.theta_eps_x,  0];

% Total misalignment transformation (body to sensor frame)
sens_eh.C_misalign = eye(3) + sens_eh.theta_eps;

% --- Non-orthogonality (detector axis non-orthogonality) ---
sens_eh.eps_xy = deg2rad(0.005);  % [rad] (~0.005 deg typical)
sens_eh.eps_xz = deg2rad(0.005);  % [rad]
sens_eh.eps_yz = deg2rad(0.005);  % [rad]

sens_eh.O = [1,              sens_eh.eps_xy, sens_eh.eps_xz;
             sens_eh.eps_xy, 1,              sens_eh.eps_yz;
             sens_eh.eps_xz, sens_eh.eps_yz, 1];

% --- Scale Factor Error ---
% Indicates error in converting measured angle to true angle
sens_eh.scale_factor_error = 0.001;  % 0.1% typical for IR sensors

% --- Temperature Sensitivity ---
sens_eh.temp_coeff = deg2rad(0.01)/10;  % [rad/°C] (0.01 deg per 10°C)
sens_eh.temp_nominal = 20;              % [°C] nominal operating temperature

% --- Sampling and Digital Conversion ---
sens_eh.update_rate = 10;                      % [Hz] typical for static sensors
sens_eh.Ts          = 1 / sens_eh.update_rate; % sampling time [s]

% ADC characteristics
sens_eh.nbits = 16;                             % 16-bit ADC (typical)
sens_eh.FS    = deg2rad(sens_eh.fov_fine);      % full scale = fine FOV [rad]
sens_eh.LSB   = sens_eh.FS / (2^sens_eh.nbits); % quantization step [rad]

% --- Earth Model Parameters (for horizon sensing) ---
h_atm = 40 * 10^3; % atmospheric height for IR horizon [m] (~40 km for 15 μm CO2)

% Effective Earth radius (including atmosphere)
R_earth_eff = R_E + h_atm;

% --- Operational Limits ---
% Minimum altitude for proper operation (FOV must see horizon)
sens_eh.h_min = 200e3;  % [m] 200 km minimum altitude

% Maximum sun angle for operation (sun interference)
sens_eh.sun_angle_max = deg2rad(30);  % [rad] 30 deg from Earth limb

% --- Noise Characteristics ---
% For Simulink Band-Limited White Noise block
sens_eh.noise_power = sens_eh.sigma_1axis^2;  % [rad^2]

% --- Random Walk (for modeling long-term drift) ---
% This is not typically specified for horizon sensors but can be estimated
sens_eh.random_walk = deg2rad(0.01) * sqrt(sens_eh.Ts);  % [rad/√s]

% Output flags
sens_eh.validity_threshold = deg2rad(80);  % Max angle from nadir for valid measurement [rad]

%% sensors - sun
T_sun = 60 * 60 * 24 * 365.25;      % solar orbit period [s]
n_sun = 2*pi / T_sun;               % solar average rotation rate [rad s^-1]
R_sun = 1.496e8;                    % solar orbit radius [km]
e_sun = deg2rad(23.45);             % solar eccentricity [-]

%% disturbances - magnetism (order 5)
j = [0.01; 0.05; 0.01;]; % magnetic dipole moment [amp m^2]
% we chose order 5 because order 6 had less than a 1% difference in the 
% total magnitude of the magnetic field vector (see mag_order_pathfinder.m)

% to have a real model of the magnetic field at each point in the orbit, we
% must simulate it beforehand.
% knowing the initial conditions of the orbit, we can simulate the rest of
% it and use magnetic field models to predict the state of the field when
% the spacecraft will be there.

% this simulation of the magnetic field is carried out in the simulation
% section at the end of this script.

%% disturbances - gravity gradient
% see simulink model

%% simulation
% simulation options
sim_step = 0.01; sim_time = round(T/15, 0);
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = num2str(sim_step);
sim_options.StartTime = "0";
sim_options.StopTime = num2str(sim_time);

% --- magnetic field model ---
% time span
n_steps = sim_time/sim_step;
tspan = linspace(0, sim_time, n_steps);

B_N = zeros(n_steps, 3);
G = zeros(5, 5 + 1);
H = zeros(5, 5 + 1);

G(5, 1) = -234.42;
G(5, 2) = 47.52;    H(5, 2) = 47.52;
G(5, 3) = 208.36;   H(5, 3) = 208.36;
G(5, 4) = -121.43;  H(5, 4) = -121.43;
G(5, 5) = 32.09;    H(5, 5) = 32.09;
G(5, 6) = 13.98;    H(5, 6) = 99.14;

% --- orbit propogation ---
% initial state vector
[r0, v0] = kep2car_SAD( ...
    a, ...
    e, ...
    inc, ...
    Omega, ...
    omega, ...
    TA0);
y0 = [r0, v0];

% set solver conditions
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% integrate
[~, Y] = ode113(@(t,y) ode_2bp_SAD(t,y,mu_E), tspan, y0, options);
r_mags = vecnorm(Y(:, 1:3), 2, 2); % create vector of position vector magnitudes

[long, lat] = ground_track_SAD(Y, T, w_E); % longitude and latitude [rad]
long = rad2deg(long); lat = rad2deg(lat); % longitude and latitude [deg]
did.disp_prog_25 = 0; did.disp_prog_50 = 0; did.disp_prog_75 = 0;

B_N = get_mag_order_5_vectorised( ...
    R_E, r_mags, ... % using r_mags instead of a because eccentricity != 0
    lat, long, ...
    G, H, ...
    5);

sim_data.time = tspan;
sim_data.signals.values = B_N;
sim_data.signals.dimensions = width(B_N);

% simulation outputs
disp("running sim")
simout = sim("project_review_simulink.slx", sim_options);
disp("sim complete")

%% references
% 1: earth-horizon dimensions
% https://satcatalog.s3.amazonaws.com/components/25/SatCatalog_-_Adcole_
% Maryland_Aerospace_-_AI-SES_IR_Earth_Sensor_-_Datasheet.pdf?lastmod=20210708041438

% 2: earth-horizon misalignment error
% https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=3096&context=smallsat

% 3: earth-horizon spec sheet
% https://electronics.leonardo.com/documents/16277707/18404907/IRES_NE_
% Attitude_Control_Sensors_LQ_mm07787_.pdf?t=1538987566453
