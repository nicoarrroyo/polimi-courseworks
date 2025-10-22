clear; close all; clc;

%% initial conditions
% constants
G = 6.67e-11;
M_e = 5.97e24;
R = 6371 + 400;
mu_E = astroConstants(13);

% inertia
ix = 0.04;
iy = 0.06;
iz = 0.08;
I = diag([ ix iy iz ]);

% angular velocities
wx0 = 0;
wy0 = 0;

n = sqrt(G*M_e/(R^3));
wz0 = n;
w0 = [wx0 wy0 wz0];

wLN = [0; 0; n;];

% DCM
ABN0 = eye(3);
ANL0 = ABN0;

c = [1; 0; 0;];

% simulation options
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.05";
sim_options.StartTime = "0";
sim_options.StopTime = "50";

%% outputs
simout = sim("lab_6_task_2_simulink", sim_options);

% time
t = simout.tout;

% angular velocity
w = simout.w.Data;
wx = w(:, 1);
wy = w(:, 2);
wz = w(:, 3);

wBL = simout.wBL.Data;

% angular acceleration
wdot = simout.wdot.Data;
wxdot = wdot(:, 1);
wzdot = wdot(:, 2);
wydot = wdot(:, 3);

% DCM - body
ABN = simout.ABN.Data;

% DCM - desired



