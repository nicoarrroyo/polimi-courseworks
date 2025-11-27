clear; close all; clc;

%% initial conditions
% constants
R = 6371 + 400; % orbit radius [km]
mu = astroConstants(13); % earth gravitational parameter [km^3 kg^-1 s^-2]

% inertia
ix = 0.04;
iy = 0.06;
iz = 0.08;
I = diag([ ix iy iz ]);

% angular velocities
wx0 = 10^-6;
%wx0 = 0;
wy0 = 10^-6;
%wy0 = 0;

n = sqrt(mu/(R^3));
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
sim_options.FixedStep = "0.1";
sim_options.StartTime = "0";
sim_options.StopTime = "20000";

%% outputs
disp("running sim")
simout = sim("lab_6_task_2_simulink", sim_options);
disp("sim complete, gathering results")

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

% DCM - LVLH
ALN = simout.ALN.Data;

% DCM - error
ABL = simout.ABL.Data;
disp("results gathered")
%animate_frame(ABN, t, "SpeedUp", str2num(sim_options.StopTime)/10)
% figure("Name", "Error DCM")




