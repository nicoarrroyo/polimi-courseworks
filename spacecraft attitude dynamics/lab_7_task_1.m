clear; close all; clc;

%% initial conditions
% constants
G = 6.67e-11;
M_e = 5.97e24;
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
R0 = (a * (1 - e^2) / (1 + (e * cos(TA0))));
n = sqrt(G*M_e/(R0^3));

% angular velocities
wx0 = 0;
wy0 = 0;
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
sim_options.FixedStep = "0.005";
sim_options.StartTime = "0";
sim_options.StopTime = "10";

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

wBL = simout.wBL.Data;

% angular acceleration
wdot = simout.wdot.Data;
wxdot = wdot(:, 1);
wzdot = wdot(:, 2);
wydot = wdot(:, 3);

% DCM - body
ABN = simout.ABN.Data;
animate_frame(ABN, t, "SpeedUp", str2double(sim_options.StopTime)/10)

% DCM - LVLH
ALN = simout.ALN.Data;

% DCM - error
ABL = simout.ABL.Data;
% figure("Name", "Error DCM")
% hold on
% title("Attitude Error DCM")
% grid on
% hold off