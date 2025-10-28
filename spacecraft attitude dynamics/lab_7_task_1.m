clear; close all; clc;

%% initial conditions
% constants
mu_E = astroConstants(13)*10^9;

% inertia
ix = 0.04;
iy = 0.06;
iz = 0.08;
I = diag([ ix iy iz ]);

% position
a = 300; % SMA [km]
e = 0; % eccentricity []
TA0 = 0; % initial true anomaly
inc = 0; % inclination [deg]
n = sqrt(mu_E/(a^3)); % average rotational rate
T = 2*pi * n; % period [s]

% sun
n_sun = 2*pi / (60*60*24*365); % average rotational rate (sun)
e_sun = 23.45; % eccentricity (sun) []

% angular velocities
w0 = [0.45; 0.52; 0.55;];

w_LN = [0; 0; n;];

% DCM
A_BN0 = eye(3);
A_NL0 = A_BN0;

c = [1; 0; 0;];

% simulation options
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.05";
sim_options.StartTime = "0";
sim_options.StopTime = num2str(T);

%% outputs
disp("running sim")
simout = sim("lab_7_task_1_simulink", sim_options);
disp("sim complete")

% time
t = simout.tout;

% angular velocity
w = simout.w.Data;
wx = w(:, 1);
wy = w(:, 2);
wz = w(:, 3);

% angular acceleration
wdot = simout.wdot.Data;
wxdot = wdot(:, 1);
wzdot = wdot(:, 2);
wydot = wdot(:, 3);

% DCM - body
A_BN = simout.A_BN.Data;
animate_frame(A_BN, t, "SpeedUp", str2double(sim_options.StopTime)/10)
