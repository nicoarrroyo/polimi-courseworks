clear
clc

% task 2
% verify the analytical solution for the symmetric case
% does the analytic solution provide a good approximation for the 
% asymmetric case?

%% set initial conditions
ix = 0.0504; % kg m^2
iy = 0.0504; % kg m^2
iz = 0.0109; % kg m^2

% defined but it could be any
wx0 = 0.1; % rad s^-1
wy0 = 0.1; % rad s^-1
wz0 = 0.1; % rad s^-1

t_end = 200;
t_step = 0.1;
t = 0:t_step:t_end;

%% analytical solution
lambda = ((iz - ix) / ix) * wz0;
wx_t = wx0*cos(lambda*t) - wy0*sin(lambda*t);
wy_t = wx0*sin(lambda*t) + wy0*cos(lambda*t);
wz_t = wz0 + 0*t;

%% get simulation outputs
simout = sim('lab_3_task_2_simulink');
wxsim = simout.wx;
wysim = simout.wy;
wzsim = simout.wz;

%% plot data
figure()
plot(wxsim, "g")
xlabel("Time (s)")
ylabel("Angular Acceleration (rad s^-2)")
title("Rotational Motion fo a 3U Cubesat")
grid on
hold on
plot(wysim, "b")
plot(wzsim, "r")

plot(t, wx_t, 'r--')
plot(t, wy_t, 'g--')
plot(t, wz_t, 'b--')
legend("wxsim", "wysim", "wzsim", "wx_t", "wy_t", "wz_t")
hold off

%% ASYMMETRIC CASE
ix_a = 0.070; % kg m^2
iy_a = 0.055; % kg m^2
iz_a = 0.025; % kg m^2

%% Define initial conditions for the asymmetric case
wx0_a = wx0; % rad s^-1
wy0_a = wy0; % rad s^-1
wz0_a = wz0; % rad s^-1

%% Calculate angular accelerations for the asymmetric case
wxdot0_a = ((iy_a - iz_a) / ix_a) * wy0_a * wz0_a;
wydot0_a = ((iz_a - ix_a) / iy_a) * wz0_a * wx0_a;
wzdot0_a = ((ix_a - iy_a) / iz_a) * wx0_a * wy0_a;

%% Calculate analytical solution for the asymmetric case
lambda_a = ((iz_a - ix_a) / ix_a) * wz0_a;
wx_t_a = wx0_a*cos(lambda_a*t) - wy0_a*sin(lambda_a*t);
wy_t_a = wx0_a*sin(lambda_a*t) + wy0_a*cos(lambda_a*t);
wz_t_a = wz0_a + 0*t;

%% Get simulation outputs for the asymmetric case
simout_a = sim('lab_3_task_2_simulink');
wxsim_a = simout_a.wx;
wysim_a = simout_a.wy;
wzsim_a = simout_a.wz;

%% plot data for the asymmetric case
figure()
plot(wxsim_a, "r")
xlabel("Time (s)")
ylabel("Angular Acceleration (rad s^-2)")
title("Rotational Motion for Asymmetric Case")
grid on
hold on
plot(wysim_a, "g")
plot(wzsim_a, "b")

plot(t, wx_t_a, 'r--')
plot(t, wy_t_a, 'g--')
plot(t, wz_t_a, 'b--')
legend("wxsim_a", "wysim_a", "wzsim_a", "wx_t_a", "wy_t_a", "wz_t_a")
hold off
