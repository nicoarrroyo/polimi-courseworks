clear
clc

% task 2
% verify the analytical solution for the symmetric case
% does the analytic solution provide a good approximation for the 
% asymmetric case?

%% symmetric case
% set initial conditions
ix = 0.0504; % kg m^2
iy = 0.0504; % kg m^2
iz = 0.0109; % kg m^2

I = diag([ix iy iz]);

w = [0.45 0.52 0.55]; % rad s^-1

% simulation options
sim_options.SolverType = "Fixed-step";
sim_options.Solver = "ode4";
sim_options.FixedStep = "0.1";
sim_options.StartTime = "0";
sim_options.StopTime = "50";

% analytical solution
wx0 = w(1); % Initial angular velocity for x-axis
wy0 = w(2); % Initial angular velocity for y-axis
wz0 = w(3); % Initial angular velocity for z-axis
t = 0:0.1:50;

lambda = ((iz - ix) / ix) * wz0;
wx_t = wx0*cos(lambda*t) - wy0*sin(lambda*t);
wy_t = wx0*sin(lambda*t) + wy0*cos(lambda*t);
wz_t = wz0 + 0*t;

% get simulation outputs
simout = sim('lab_3_task_2_simulink', sim_options);
t = simout.tout;
wsim = simout.w.Data;
wxsim = wsim(:, 1);
wysim = wsim(:, 2);
wzsim = wsim(:, 3);

% plot data
figure()
plot(t, wxsim, "g")
xlabel("Time (s)")
ylabel("Angular Acceleration (rad s^-2)")
title("Rotational Motion fo a 3U Cubesat")
grid on
hold on
plot(t, wysim, "b")
plot(t, wzsim, "r")

plot(t, wx_t, 'r--')
plot(t, wy_t, 'g--')
plot(t, wz_t, 'b--')
legend("wxsim", "wysim", "wzsim", "wx_t", "wy_t", "wz_t")
hold off

%% ASYMMETRIC CASE
ix_a = 0.070; % kg m^2
iy_a = 0.055; % kg m^2
iz_a = 0.025; % kg m^2

I_a = diag([ix_a iy_a iz_a]);

% Define initial conditions for the asymmetric case
w_a = [0.45 0.52 0.55]; % rad s^-1

% Calculate analytical solution for the asymmetric case
wx0_a = w_a(1); % Initial angular velocity for x-axis
wy0_a = w_a(2); % Initial angular velocity for y-axis
wz0_a = w_a(3); % Initial angular velocity for z-axis

lambda_a = ((iz_a - ix_a) / ix_a) * wz0_a;
wx_t_a = wx0_a*cos(lambda_a*t) - wy0_a*sin(lambda_a*t);
wy_t_a = wx0_a*sin(lambda_a*t) + wy0_a*cos(lambda_a*t);
wz_t_a = wz0_a + 0*t;
w_t_a = [wx_t_a wy_t_a wz_t_a];

% Get simulation outputs for the asymmetric case
simout_a = sim('lab_3_task_2_a_simulink', sim_options);
wsim_a = simout_a.w_a.Data;
t_a = simout.tout;
wxsim_a = wsim_a(:, 1);
wysim_a = wsim_a(:, 2);
wzsim_a = wsim_a(:, 3);

% plot data for the asymmetric case
figure()
plot(t, wxsim_a, "r")
xlabel("Time (s)")
ylabel("Angular Acceleration (rad s^-2)")
title("Rotational Motion for Asymmetric Case")
grid on
hold on
plot(t, wysim_a, "g")
plot(t, wzsim_a, "b")

plot(t, wx_t_a, 'r--')
plot(t, wy_t_a, 'g--')
plot(t, wz_t_a, 'b--')
legend("wxsim_a", "wysim_a", "wzsim_a", "wx_t_a", "wy_t_a", "wz_t_a")
hold off
