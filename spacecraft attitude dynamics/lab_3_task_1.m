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

%% for verification
wxdot0 = ((iy - iz) / ix) * wy0 * wz0;
wydot0 = ((iz - ix) / iy) * wz0 * wx0;
wzdot0 = ((ix - iy) / iz) * wx0 * wy0;

%% get simulation outputs
simout = sim('lab_3_task_1_simulink');
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
