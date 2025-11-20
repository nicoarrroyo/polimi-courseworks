clear
clc
close all

% task 4
% test the stability of a dual-spin spacecraft
% evaluate stability by playing with different values of the off-axis 
% initial angular velocities and wheel velocity
% do the same with iz = 0.0504, iy = 0.0109. reason the results

%% case 1
% set initial conditions
ix = 0.070; % kg m^2
iy = 0.055; % kg m^2
iz = 0.025; % kg m^2
ir = 0.005; % kg m^2

wx0 = 1e-6; % rad s^-1
wy0 = 1e-6; % rad s^-1
wz0 = 0.02; % rad s^-1
wr0 = 2*pi; % rad s^-1

% simulation options
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.1";
sim_options.StartTime = "0";
sim_options.StopTime = "20";

% get simulation outputs
simout = sim('lab_3_task_4_simulink', sim_options);
time = simout.tout;
wxdot = simout.wxdot.Data;
wydot = simout.wydot.Data;
wzdot = simout.wzdot.Data;
wx = simout.wx.Data;
wy = simout.wy.Data;
wz = simout.wz.Data;

% plot data
figure()
plot(time, wxdot, "r")
xlabel("Time (s)")
ylabel("Angular Acceleration (rad s^-^2)")
title("Rotational Motion for a 3U Cubesat (ang. accel.)")
grid on
hold on
plot(time, wydot, "g")
plot(time, wzdot, "b")
legend("wxdot", "wydot", "wzdot")
hold off

figure()
plot(time, wx, "r")
xlabel("Time (s)")
ylabel("Angular Velocity (rad s^-^1)")
title("Rotational Motion for a 3U Cubesat (ang. vel.)")
grid on
hold on
plot(time, wy, "g")
plot(time, wz, "b")
legend("wx", "wy", "wz")
hold off


%% case 2
% note: simulink will keep the latest values stored in the worksspace
% set initial conditions
iz = 0.0504; % kg m^2
iy = 0.0109; % kg m^2

% simulation options
% smaller step size needed to avoid singularity
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.01";
sim_options.StartTime = "0";
sim_options.StopTime = "20";

% get simulation outputs
simout = sim('lab_3_task_4_simulink', sim_options);
time = simout.tout;
wxdot = simout.wxdot.Data;
wydot = simout.wydot.Data;
wzdot = simout.wzdot.Data;
wx = simout.wx.Data;
wy = simout.wy.Data;
wz = simout.wz.Data;

% plot data
figure()
plot(time, wxdot, "r")
xlabel("Time (s)")
ylabel("Angular Acceleration (rad s^-^2)")
title("Rotational Motion for a 3U Cubesat (ang. accel.)")
grid on
hold on
plot(time, wydot, "g")
plot(time, wzdot, "b")
legend("wxdot", "wydot", "wzdot")
hold off

figure()
plot(time, wx, "r")
xlabel("Time (s)")
ylabel("Angular Velocity (rad s^-^1)")
title("Rotational Motion for a 3U Cubesat (ang. vel.)")
grid on
hold on
plot(time, wy, "g")
plot(time, wz, "b")
legend("wx", "wy", "wz")
hold off
