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
Omega = 0;                        % Right Ascension of Ascending node [rad]
omega = 0;
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

%% Sensors - Gyroscope
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

% DCM GYRO BODY (fixed matrix, doesnt change in time, in our case it is the
% identity matrix)
sens_gyro_dcm_body = eye(3);

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
sens_gyro_LSB = sens_gyro_FS / (2 ^ sens_gyro_nbits); % Least significant bit

%% Sensors - Magnetometer
% Based on AAC SpaceQuest MAG-3 3-axis fluxgate magnetometer
% Space-qualified, TRL-9, flight-proven on numerous missions

% --- Physical Dimensions ---
sens_mag_length = 82.6e-3;  % length [m]
sens_mag_width  = 35.1e-3;  % width [m]
sens_mag_height = 32.3e-3;  % height [m]
sens_mag_m      = 100e-3;   % mass [kg]

% --- Performance Specifications ---
% Accuracy (as percentage of full scale)
sens_mag_accuracy = 0.0075;  % ±0.75% of FS (0.5% typical)

% Linearity
sens_mag_linearity = 0.00015;  % ±0.015% of FS

% Sensitivity (output voltage per field strength)
sens_mag_sensitivity = 100e-6;  % 100 μV/nT [V/nT]

% Field measurement range
sens_mag_range = 65e-6;  % ±65 μT (standard range) [T]

% Noise density
sens_mag_noise_density = 12e-12;  % 12 pT/√Hz @1Hz [T/√Hz]

% Noise density (alternative for 5V model)
sens_mag_noise_density_5v = 100e-12;  % <100 pT/√Hz @1Hz [T/√Hz]

% Zero field offset
sens_mag_zero_field_offset = 25e-3;  % ±0.025 V analog output

% --- Temperature Effects ---
% Scale factor temperature shift
sens_mag_scale_factor_temp = 0.00007;  % 0.007% FS/°C

% Zero shift with temperature
sens_mag_zero_temp_shift = 0.6e-9;  % ±0.6 nT/°C [T/°C]

% Operating temperature range
sens_mag_temp_min = -55;  % [°C]
sens_mag_temp_max = 85;   % [°C]
sens_mag_temp_nominal = 20;  % [°C]

% --- Axial Alignment ---
% Orthogonality between axes
sens_mag_orthogonality_error = deg2rad(1);  % Better than ±1 deg [rad]

% --- Perming Effects (magnetic memory) ---
% Susceptibility to perming (residual magnetization)
sens_mag_perming_shift = 8e-9;  % ±8 nT shift with ±5 Gauss applied [T]
sens_mag_perming_field = 5e-4;  % ±5 Gauss = 500 μT [T]

% --- Frequency Response ---
sens_mag_bandwidth = 500;  % 3dB @ >500 Hz (up to 4 kHz wideband)

% --- Sampling Rate ---
sens_mag_update_rate = 10;  % [Hz] typical for attitude determination
sens_mag_Ts = 1 / sens_mag_update_rate;  % sampling time [s]

% --- Electrical Specifications ---
sens_mag_voltage_input = [15, 34];  % 15-34 VDC (or 5V regulated)
sens_mag_current = 30e-3;  % 30 mA at any input voltage [A]
sens_mag_power = sens_mag_voltage_input(1) * sens_mag_current;  % ~0.45 W

% Analog output options
sens_mag_output_voltage = 10;  % ±10 V output [V]
sens_mag_output_range_T = 100e-6;  % corresponds to ±100 μT [T]

% Alternative: ±5V output for ±60μT
% sens_mag_output_voltage = 5;
% sens_mag_output_range_T = 60e-6;

% --- Digital Conversion (ADC) ---
sens_mag_nbits = 16;  % 16-bit ADC typical for spacecraft magnetometers
sens_mag_LSB = (2 * sens_mag_range) / (2^sens_mag_nbits);  % [T]

% --- Derived Specifications ---
% RMS noise at update rate
sens_mag_noise_rms = sens_mag_noise_density * sqrt(sens_mag_update_rate);  % [T]

% 1-sigma measurement uncertainty (from accuracy spec)
sens_mag_sigma_1axis = sens_mag_accuracy * sens_mag_range;  % [T]

% Scale factor (converts measured field to true field)
sens_mag_scale_factor = 1.0;  % nominal, will add error term

% --- Misalignment Matrix ---
% Installation misalignment (sensor to body frame)
sens_mag_theta_eps_x = deg2rad(0.5);  % [rad] typical mounting tolerance
sens_mag_theta_eps_y = deg2rad(0.5);  % [rad]
sens_mag_theta_eps_z = deg2rad(0.5);  % [rad]

sens_mag_theta_eps = [0,                     -sens_mag_theta_eps_z,  sens_mag_theta_eps_y;
                      sens_mag_theta_eps_z,   0,                    -sens_mag_theta_eps_x;
                     -sens_mag_theta_eps_y,   sens_mag_theta_eps_x,  0];

% Total misalignment transformation
sens_mag_C_misalign = eye(3) + sens_mag_theta_eps;

% --- Non-Orthogonality Matrix (between sensor axes) ---
% Based on ±1° orthogonality spec
sens_mag_eps_xy = deg2rad(0.5);  % [rad]
sens_mag_eps_xz = deg2rad(0.5);  % [rad]
sens_mag_eps_yz = deg2rad(0.5);  % [rad]

sens_mag_O = [1,                 sens_mag_eps_xy, sens_mag_eps_xz;
              sens_mag_eps_xy,   1,               sens_mag_eps_yz;
              sens_mag_eps_xz,   sens_mag_eps_yz, 1];

% % --- Hard Iron and Soft Iron Effects ---
% % Hard iron: constant magnetic field from spacecraft (ferromagnetic materials)
% sens_mag_hard_iron = [50e-9; 30e-9; -40e-9];  % [T] typical spacecraft disturbance
% 
% % Soft iron: induced magnetization (scaling and rotation)
% % Represented as a 3x3 matrix deviation from identity
% sens_mag_soft_iron = eye(3) + 0.01 * randn(3);  % 1% soft iron effect

% --- Bias Instability (Random Walk) ---
% Not typically specified for magnetometers, but exists
sens_mag_bias_instability = 0.1e-9;  % [T/s] estimated (0.1 nT/s)

% --- Earth's Magnetic Field (for reference) ---
% These values are location-dependent; typical LEO values
B_earth_magnitude_LEO = 30e-6;  % ~30 μT typical in LEO [T]
B_earth_inclination = deg2rad(60);  % typical inclination [rad]

% --- Output Flags ---
sens_mag_valid_threshold = 1.5 * sens_mag_range;  % saturation threshold

% --- Radiation Tolerance ---
sens_mag_radiation_TID = 10e3;  % >10 krad Total Ionizing Dose [rad]

%% Sensors - Earth Horizon
% Based on MAI-SES IR Earth Sensor specifications

% --- Physical Dimensions ---
sens_eh_m      = 33e-3;     % mass [kg]
sens_eh_length = 43.3e-3;   % length [m] - CORRECTED to meters
sens_eh_width  = 31.8e-3;   % width [m]
sens_eh_height = 20.7e-3;   % height [m]

% --- Optical Characteristics ---
sens_eh_fov_coarse       = 60;      % coarse field of view [deg]
sens_eh_fov_fine         = 7;       % fine field of view [deg]
sens_eh_fov_coarse_rad   = deg2rad(sens_eh_fov_coarse);  % [rad]
sens_eh_fov_fine_rad     = deg2rad(sens_eh_fov_fine);    % [rad]

% --- Accuracy Specifications (using fine FOV resolution) ---
sens_eh_resolution_fine   = 0.25;    % fine FOV resolution [deg]
sens_eh_resolution_coarse = 1.0;     % coarse FOV resolution [deg]

% Convert to radians for simulation
sens_eh_sigma_1axis = deg2rad(sens_eh_resolution_fine); % 1-sigma noise per axis [rad]

% For 2-axis measurement (pitch and roll from nadir)
sens_eh_sigma_pitch = sens_eh_sigma_1axis;  % [rad]
sens_eh_sigma_roll  = sens_eh_sigma_1axis;  % [rad]

% --- Bias and Systematic Errors ---
% Static bias (offset from true nadir direction)
sens_eh_bias_pitch = deg2rad(0.1);   % pitch bias [rad] (~0.1 deg typical)
sens_eh_bias_roll  = deg2rad(0.1);   % roll bias [rad]

% Bias instability (slow drift over time)
sens_eh_bias_instability = deg2rad(0.05/3600); % [rad/s] (0.05 deg/hr typical)

% --- Misalignment Errors ---
% Installation misalignment (sensor to body frame)
sens_eh_theta_eps_x = deg2rad(0.01);  % [rad] (~0.01 deg typical mounting tolerance)
sens_eh_theta_eps_y = deg2rad(0.01);  % [rad]
sens_eh_theta_eps_z = deg2rad(0.01);  % [rad]

% Misalignment matrix (small angle approximation)
% this is applicable because we won't use earth-horizon for detumbling
sens_eh_theta_eps = [0,                    -sens_eh_theta_eps_z,  sens_eh_theta_eps_y;
                     sens_eh_theta_eps_z,   0,                   -sens_eh_theta_eps_x;
                    -sens_eh_theta_eps_y,   sens_eh_theta_eps_x,  0];

% Total misalignment transformation (body to sensor frame)
sens_eh_C_misalign = eye(3) + sens_eh_theta_eps;

% --- Non-orthogonality (detector axis non-orthogonality) ---
sens_eh_eps_xy = deg2rad(0.005);  % [rad] (~0.005 deg typical)
sens_eh_eps_xz = deg2rad(0.005);  % [rad]
sens_eh_eps_yz = deg2rad(0.005);  % [rad]

sens_eh_O = [1,              sens_eh_eps_xy, sens_eh_eps_xz;
             sens_eh_eps_xy, 1,              sens_eh_eps_yz;
             sens_eh_eps_xz, sens_eh_eps_yz, 1];

% --- Scale Factor Error ---
% Indicates error in converting measured angle to true angle
sens_eh_scale_factor_error = 0.001;  % 0.1% typical for IR sensors

% --- Temperature Sensitivity ---
sens_eh_temp_coeff = deg2rad(0.01)/10;  % [rad/°C] (0.01 deg per 10°C)
sens_eh_temp_nominal = 20;              % [°C] nominal operating temperature

% --- Sampling and Digital Conversion ---
sens_eh_update_rate = 10;                      % [Hz] typical for static sensors
sens_eh_Ts          = 1 / sens_eh_update_rate; % sampling time [s]

% ADC characteristics
sens_eh_nbits = 16;                            % 16-bit ADC (typical)
sens_eh_FS    = deg2rad(sens_eh_fov_fine);     % full scale = fine FOV [rad]
sens_eh_LSB   = sens_eh_FS / (2^sens_eh_nbits);% quantization step [rad]

% --- Earth Model Parameters (for horizon sensing) ---
R_earth = 6371e3;  % Earth radius [m] (mean)
h_atm   = 40e3;    % atmospheric height for IR horizon [m] (~40 km for 15 μm CO2)

% Effective Earth radius (including atmosphere)
R_earth_eff = R_earth + h_atm;

% --- Operational Limits ---
% Minimum altitude for proper operation (FOV must see horizon)
sens_eh_h_min = 200e3;  % [m] 200 km minimum altitude

% Maximum sun angle for operation (sun interference)
sens_eh_sun_angle_max = deg2rad(30);  % [rad] 30 deg from Earth limb

% --- Noise Characteristics ---
% White noise power spectral density
sens_eh_noise_psd = sens_eh_sigma_1axis / sqrt(sens_eh_Ts);  % [rad/√Hz]

% For Simulink Band-Limited White Noise block
sens_eh_noise_power = sens_eh_sigma_1axis^2;  % [rad^2]

% --- Random Walk (for modeling long-term drift) ---
% This is not typically specified for horizon sensors but can be estimated
sens_eh_random_walk = deg2rad(0.01) * sqrt(sens_eh_Ts);  % [rad/√s]

% Output flags
sens_eh_validity_threshold = deg2rad(80);  % Max angle from nadir for valid measurement [rad]

%% sensors - sun
T_sun = 60 * 60 * 24 * 365.25;      % solar orbit period [s]
n_sun = 2*pi / T_sun;               % solar average rotation rate [rad s^-1]
R_sun = 1.496e8;                    % solar orbit radius [km]
e_sun = deg2rad(23.45);             % solar eccentricity [-]

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
sim_step = 0.01; sim_time = T;
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = num2str(sim_step);
sim_options.StartTime = "0";
sim_options.StopTime = num2str(round(sim_time, 0));

% magnetic field model
% time span
period = 2*pi * sqrt(a^3 / mu_E); % orbital period [s]
tspan = linspace(0, sim_step, sim_time);

% --- orbit propogation ---
% initial state vector
[r0, v0] = kep2car( ...
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
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y0, options);

[long, lat] = ground_track(Y, T, w_E); % longitude and latitude [radians]

for i = 1:(sim_time/sim_step)
    
end
B_N = get_mag_order_5(R_E, a, rad2deg(lat), rad2deg(long));

% simulation outputs
disp("running sim")
%simout = sim("project_review_simulink.slx", sim_options);
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
