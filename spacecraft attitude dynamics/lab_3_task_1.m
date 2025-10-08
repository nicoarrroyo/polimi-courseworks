clear
clc

% task 1
% simulate the rotational motion of a 3u cubesat
% use a scope to analyze the output
% plot the output in MATLAB with proper axis labels and units

%% set initial conditions
ix = 0.070; % kg m^2
iy = 0.055; % kg m^2
iz = 0.025; % kg m^2

wx0 = 0.45; % rad s^-1
wy0 = 0.52; % rad s^-1
wz0 = 0.55; % rad s^-1

%% simulation options
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.1";
sim_options.StartTime = "0";
sim_options.StopTime = "100";

%% get simulation outputs
simout = sim('lab_3_task_1_simulink', sim_options);
wx = simout.wx;
wy = simout.wy;
wz = simout.wz;

%% plot data
figure()
plot(wx)
xlabel("Time (s)")
ylabel("Angular Velocity (rad s^-1)")
title("Rotational Motion for a 3U Cubesat")
grid on
hold on
plot(wy)
plot(wz)
legend("wx", "wy", "wz")
hold off
