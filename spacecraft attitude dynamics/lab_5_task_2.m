clear
%clc
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
w0 = [wx0; wy0; wz0;];
h_body = I * w0;

h_ref = [1; 0; 0;];

% A312 = @(phi, theta, psi) [
%         cos(psi)*cos(phi) + sin(psi)*sin(phi)*sin(theta), ...
%         -cos(psi)*sin(phi) + sin(psi)*cos(phi)*sin(theta), ...
%         sin(psi)*cos(theta);
%         sin(phi)*cos(theta), ...
%         cos(phi)*cos(theta), ...
%         -sin(theta);
%         -sin(psi)*cos(phi) + cos(psi)*sin(phi)*sin(theta), ...
%         sin(psi)*sin(phi) + cos(psi)*cos(phi)*sin(theta), ...
%         cos(theta)*cos(psi)];

A312 = @(phi, theta, psi) [ % all signs changed
        cos(psi)*cos(phi) - sin(psi)*sin(phi)*sin(theta), ...
        cos(psi)*sin(phi) + sin(psi)*cos(phi)*sin(theta), ...
        -sin(psi)*cos(theta);
        -sin(phi)*cos(theta), ...
        cos(phi)*cos(theta), ...
        sin(theta);
        sin(psi)*cos(phi) + cos(psi)*sin(phi)*sin(theta), ...
        sin(psi)*sin(phi) - cos(psi)*cos(phi)*sin(theta), ...
        cos(theta)*cos(psi);];

fun = @(euler_angles) h_ref - A312(euler_angles(1), euler_angles(2), euler_angles(3))' * h_body;
guess = [0; 0; 0];
options = optimoptions("fsolve", "Display", "iter-detailed", ...
    "TolFun", 1e-10, "TolX", 1e-10);
euler_angles0 = fsolve(fun, guess, options);

phi0 = euler_angles0(1);
theta0 = euler_angles0(2);
psi0 = euler_angles0(3);
A0 = A312(phi0, theta0, psi0);

h_N0 = A0' * h_body;

h_error0 = cross(h_ref, h_N0);
disp("angular momentum error")
disp(h_error0)

%% simulation options
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.1";
sim_options.StartTime = "0";
sim_options.StopTime = "15";

%% run sim
simout = sim("lab_5_task_2_simulink", sim_options);
t = simout.tout;
w = simout.w.Data;
wx = w(:, 1);
wy = w(:, 2);
wz = w(:, 3);

wdot = simout.wdot.Data;

phi_sim = simout.phi.Data;
psi_sim = simout.psi.Data;
theta_sim = simout.theta.Data;
angles_sim = [phi_sim, psi_sim, theta_sim];

phidot = simout.phidot.Data;
psidot = simout.psidot.Data;
thetadot = simout.thetadot.Data;
anglesdot = [phidot, psidot, thetadot];

A = simout.A.Data;

%% simulate reference frames
animate_frame(A, t);

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
verify = A(:, :, end) * A(:, :, end)';
disp("verify")
disp(verify)

%% plot angles data
figure()
plot(t, angles_sim)
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