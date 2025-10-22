clear
clc
close all

% lab 5 task 1
% Simulate the kinematics using Euler Angle
% Simulate the kinematics by accounting for the presence of singularities
% build the DCM from the Euler Angles

%% set initial conditions
ix = 0.070;
iy = 0.055;
iz = 0.025;
I = diag([ix iy iz]);

wx0 = 0.45;
wy0 = 0.52;
wz0 = 0.55;
w0 = [wx0 wy0 wz0];

A0 = eye(3);

phi0 = 0.1;
psi0 = 0.2;
theta0 = 0.3;

%% simulation options
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.01";
sim_options.StartTime = "0";
sim_options.StopTime = "100";

%% run sim
simout = sim("lab_5_task_1_simulink", sim_options);
t = simout.tout;
w = simout.w.Data;
wx = w(:, 1);
wy = w(:, 2);
wz = w(:, 3);

A = simout.A.data;

wdot = simout.wdot.Data;

phi = simout.phi.Data;
psi = simout.psi.Data;
theta = simout.theta.Data;
angles = [phi, psi, theta];

phidot = simout.phidot.Data;
psidot = simout.psidot.Data;
thetadot = simout.thetadot.Data;
anglesdot = [phidot, psidot, thetadot];

threshold = 0.5;
total_counter = 0;

%% plot w data
figure()
plot(t, wx, "r")
xlabel("Time (s)")
ylabel("Angular Velocity (rad s^-1)")
title("Rotational Motion for a 3U Cubesat")
grid on
hold on
plot(t, wy, "g")
plot(t, wz, "b")
legend("wx", "wy", "wz")
hold off

%% plot DCM A data
% figure()
% plot(t, A(3, 2, 1))
% xlabel("Time (s)")
% ylabel("Attitude Parameter A")
% title("Rotational Motion for a 3U Cubesat")
% grid on
% hold off

%% plot wdot data
figure()
plot(t, wdot)
xlabel("Time (s)")
ylabel("Angular Acceleration (rad s^-2)")
title("Rotational Motion for a 3U Cubesat")
grid on
legend("wdotx", "wdoty", "wdotz")
hold off

%% verify (task 2)
verify = A(:, :, end) * A(:, :, end)'

%% plot angles data
figure()
plot(t, angles)
xlabel("Time (s)")
ylabel("Angular Velocity (rad s^-1)")
title("Rotational Motion for a 3U Cubesat")
grid on
legend("phi", "psi", "theta")
hold off

%% plot angles dot data
figure()
plot(t, anglesdot)
xlabel("Time (s)")
ylabel("Angular Velocity (rad s^-1)")
title("Kinematics for a 3U Cubesat")
grid on
legend("phidot", "psidot", "thetadot")
hold off
