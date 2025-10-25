clear
clc
close all

% lab 5 task 2 use quiver3
%% set initial conditions
ix = 0.070;
iy = 0.055;
iz = 0.025;
I = diag([ix iy iz]);

wx0 = 0.2 + (2*pi - 0.2)*rand;
wy0 = 0.1;
wz0 = 0.1;
w0 = [wx0 wy0 wz0];

A0 = eye(3);

phi0 = 0;
psi0 = 0;
theta0 = 0;

%% simulation options
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.01";
sim_options.StartTime = "0";
sim_options.StopTime = "5";

%% run sim
simout = sim("lab_5_task_2_simulink", sim_options);
t = simout.tout;
w = simout.w.Data;
wx = w(:, 1);
wy = w(:, 2);
wz = w(:, 3);

wdot = simout.wdot.Data;

phi = simout.phi.Data;
psi = simout.psi.Data;
theta = simout.theta.Data;
angles = [phi, psi, theta];

phidot = simout.phidot.Data;
psidot = simout.psidot.Data;
thetadot = simout.thetadot.Data;
anglesdot = [phidot, psidot, thetadot];

A = simout.A.Data;

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
ylabel("Euler Angles [rad]")
title("Euler Angles (phi, theta, psi)")
grid on
legend("phi", "psi", "theta")
hold off

%% plot angles dot data
figure()
plot(t, anglesdot)
xlabel("Time (s)")
ylabel("Time Derivative of Euler Angles [rad/s]")
title("Time Derivative of Euler Angles (phi, theta, psi)")
grid on
legend("phidot", "psidot", "thetadot")
hold off
