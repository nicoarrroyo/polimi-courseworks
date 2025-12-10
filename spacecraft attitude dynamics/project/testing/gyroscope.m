clear; close all; clc;

%% Data
% manca: ortogonality matrix O , capire brown e pink noise
% nota bene: la matrice theta non si costruisce dalla [w]^, sono 
% semplicemente degli angoli iniziali scelti da noi.
% la w è una grandezza continua, la latenza lo shifta di un certo tempo
% verso destra perchè il giroscopio ci mette una certa quantità di tempo a
% processare l'informazione.
m=0.052;                  % [kg] mass of the gyroscope
length=44.8;              % [mm] 
height=38.6;              % [mm]
deep=21.5;                % [mm]

mis_err=1/1000;           % [mrad] misalignment error
run_run_bias=4/3600;      % [degrees/second] bias due to sensor turn on
static_temp_bias=9/3600;  % [degrees/second] bias due to static temperature
SFE = 500*1e-6;           % [-] scale factor
SFN = 15*1e-6;            % [-] non linearity scale factor
update_rate = 100;        % [Hz]
Ts= 1/update_rate;        % [s]  sampling time
D_bias_inst=0.3/3600;     % [degrees/second] bias instability
D_ARW=0.15/sqrt(3600);    % [degrees/sqrt(second)] angle random walk
RRW=1e-3;                 % rate random walk guess
D_RW=RRW/sqrt(3600);      %  [degrees/second]
D_wn = D_ARW/sqrt(3600);  % Amplitude of white noise 
c_time=200;               % [s]
FS=10;                    % [V] full scale
nbits=24;                 % [-] # bits
LSB=FS/(2*exp(nbits));    % Least significant bit

% Keplerian parameters
a = 7000;                         % [km]
e = 0.001;                        % [-]
i = deg2rad(20);                  % [rad]
theta0 = 0;                       % [rad]

M_e = 5.97e24;                    % [kg]
G = 6.67e-20;                     % [km^3/(kg s^2)]
mu_e = G*M_e;                     % [km^3/s^2]
R_e = 6378;                       % [km]
omega_e = 7.2916e-5;              % [rad/s]
n = sqrt(mu_e/(a^3));             % [rad/s]
T = 2*pi / n;                     % [s]

% Inertia matrix
Ix = 25;                          % [kg m^2]
Iy = 35;                          % [kg m^2]
Iz = 50;                          % [kg m^2]
I = diag([Ix Iy Iz]);             % [kg m^2]


% Initial Conditions
w0 = [1e-6; 1e-6; n];             % [rad/s]
A_bn0 = eye(3);                   % [rad]
A_ln0 = A_bn0;                    % [rad]
c0 = [1; 0; 0];                   % [rad]

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


%% Simulation Options
sim_options.SolverType = 'Fixed-step';      % Set the solver type to Fixed-step
sim_options.Solver = 'ode4';                % Select ode4 as solver
sim_options.FixedStep = '0.01';              % Select a time step of 0.1 s
sim_options.StartTime = '0';                % Start from 0 seconds [default]
sim_options.StopTime = 'T';                % End the simulation at t=10s

deltat=0.01;
%% Simulate the model 

result = sim('gyroscope_sym.slx', sim_options);
