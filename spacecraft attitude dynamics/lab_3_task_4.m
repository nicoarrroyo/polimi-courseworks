clear
clc
close all

% task 4
% test the stability of a dual-spin spacecraft
% evaluate stability by playing with different values of the off-axis 
% initial angular velocities and wheel velocity
% do the same with iz = 0.0504, iy = 0.0109. reason the results

%% set initial conditions
I = diag([0.070 0.055 0.025]); % kg m^2
Ir = 0.005; % kg m^2

w0 = [1e-6 1e-6 0.02]; % rad s^-1
wr = 2*pi; % rad s^-1

sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.1";
sim_options.StartTime = "0";
sim_options.StopTime = "1000";

%% get simulation outputs
simout = sim('lab_3_task_4_simulink', sim_options);
time = simout.tout;
w = simout.w.Data;
wdot = simout.wdot.Data;

%% plot data
wxdot = wdot(:, 1);
wydot = wdot(:, 2);
wzdot = wdot(:, 3);

figure()
plot(time, wxdot)
xlabel("Time (s)")
ylabel("Angular Acceleration (rad s^-2)")
title("Rotational Motion fo a 3U Cubesat")
grid on
hold on
plot(time, wydot)
plot(time, wzdot)
legend("wxdot", "wydot", "wzdot")
hold off
