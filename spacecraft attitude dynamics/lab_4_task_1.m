clear
clc
close all

% lab 4 task 1
% simulate kinematics using the DCM with standard integration
% Simulate the kinematics starting trying different values of the inertia 
% tensor and initial conditions
% Test the orthonormality of the propagated solution

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

wsk0 = [
    [0 -wz0 wy0]
    [wz0 0 -wx0]
    [-wy0 wx0 0]];


%% simulation options
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode5";
sim_options.FixedStep = "0.1";
sim_options.StartTime = "0";
sim_options.StopTime = "100";

%% run sim
simout = sim("lab_4_task_1_simulink", sim_options);
t = simout.tout;
w = simout.w.Data;
wx = w(:, 1);
wy = w(:, 2);
wz = w(:, 3);

A = simout.A.data;

wdot = simout.wdot.Data;

%% plot data
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

%% plot data
figure()
plot(t, A(3, 2, 1))
xlabel("Time (s)")
ylabel("Attitude Parameter A")
title("Rotational Motion for a 3U Cubesat")
grid on
hold off

%% plot data
figure()
plot(t, wdot)
xlabel("Time (s)")
ylabel("Angular Acceleration (rad s^-2)")
title("Rotational Motion for a 3U Cubesat")
grid on
legend("wdotx", "wdoty", "wdotz")
hold off
