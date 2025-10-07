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
t = linspace(0, 200, 0.1);

%% for verification
wxdot0 = ((iy - iz) / ix) * wy0 * wz0;
wydot0 = ((iz - ix) / iy) * wz0 * wx0;
wzdot0 = ((ix - iy) / iz) * wx0 * wy0;

%% analytical solution
lambda = ((iz - ix) / ix) * wz0;
wx_t = wx0.*cos(lambda.*t) - wy0.*sin(lambda.*t);
wy_t = wx0*sin(lambda*t) + wy0*cos(lambda*t);
wz_t = wz0;

%% get simulation outputs
simout = sim('lab_3_task_2_simulink');
wxdot = simout.wxdot;
wydot = simout.wydot;
wzdot = simout.wzdot;

%% plot data
figure()
plot(wxdot)
xlabel("Time (s)")
ylabel("Angular Acceleration (rad s^-2)")
title("Rotational Motion fo a 3U Cubesat")
grid on
hold on
plot(wydot)
plot(wzdot)
legend("wxdot", "wydot", "wzdot")
hold off

%% verify initial conditions
disp("=== CONTROL INTIAL CONDITIONS ===")
fprintf("initial wxdot: CALC %d SIM: %d\n", wxdot0, wxdot.Data(1))
fprintf("initial wydot: CALC %d SIM: %d\n", wydot0, wydot.Data(1))
fprintf("initial wzdot: CALC %d SIM: %d\n", wzdot0, wzdot.Data(1))
disp("=== CONTROL INTIAL CONDITIONS ===")

% %% asymmetric case
% ix = 0.070; % kg m^2
% iy = 0.055; % kg m^2
% iz = 0.025; % kg m^2
