clear; close all; clc;

%% initial conditions
% constants
mu_E = astroConstants(13);

% inertia
ix = 0.070;
iy = 0.055;
iz = 0.025;
I = diag([ ix iy iz ]);

% position
a = 7093; % SMA [km]
e = 0.5; % eccentricity []
TA0 = deg2rad(0); % initial true anomaly [rad]
inc = deg2rad(98.27); % inclination [rad]
n = sqrt(mu_E/(a^3)); % average rotational rate
T = 2*pi * n; % period [s]

% sun
n_sun = 2*pi / (60*60*24*365); % average rotational rate (sun)
r_sun = 1.496e8;
e_sun = deg2rad(23.45); % eccentricity (sun) []

% angular velocities
w0 = [0.45; 0.52; 0.55;];

% DCM
A_BN0 = eye(3);

% simulation options
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.01";
sim_options.StartTime = "0";
sim_options.StopTime = "100";

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

% angular acceleration
wdot = simout.wdot.Data;
wxdot = wdot(:, 1);
wzdot = wdot(:, 2);
wydot = wdot(:, 3);

% DCM - body
A_BN = simout.A_BN.Data;

%% Plots
TA = simout.TA.Data;
r = simout.r.Data;
r_N = simout.r_N.Data;   
r_B = simout.r_B.Data;
S_N = simout.S_N.Data;   
S_B = simout.S_B.Data;

% True anomaly
figure();
plot(t, rad2deg(TA), 'LineWidth', 2);
xlabel('t [s]');
ylabel('$\theta$ [deg]', 'Interpreter', 'latex');
title('True Anomaly');
grid on;
xlim([t(1), t(end)]);

% Orbital radius
figure();
plot(t, r, 'LineWidth', 2);
xlabel('t [s]');
ylabel('r [km]');
title('Orbital Radius');
grid on;
xlim([t(1), t(end)]);
yline(a, '--r', 'Semi-major axis');

% Position in inertial frame
figure();
subplot(3,1,1);
plot(t, r_N(:,1), 'LineWidth', 2);
ylabel('$r_N^x$ [km]', 'Interpreter', 'latex');
grid on;
title('Position in Inertial Frame [km]');
subplot(3,1,2);
plot(t, r_N(:,2), 'LineWidth', 2);
ylabel('$r_N^y$ [km]', 'Interpreter', 'latex');
grid on;
subplot(3,1,3);
plot(t, r_N(:,3), 'LineWidth', 2);
ylabel('$r_N^z$ [km]', 'Interpreter', 'latex');
xlabel('t [s]');
grid on;

% Position in body frame
figure();
subplot(3,1,1);
plot(t, r_B(:,1), 'LineWidth', 2);
ylabel('$r_B^x$ [km]', 'Interpreter', 'latex');
grid on;
title('Position in Body Frame [km]');
subplot(3,1,2);
plot(t, r_B(:,2), 'LineWidth', 2);
ylabel('$r_B^y$ [km]', 'Interpreter', 'latex');
grid on;
subplot(3,1,3);
plot(t, r_B(:,3), 'LineWidth', 2);
ylabel('$r_B^z$ [km]', 'Interpreter', 'latex');
xlabel('t [s]');
grid on;
