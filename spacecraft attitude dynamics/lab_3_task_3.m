clear
clc

% task 3
% Assess the stability of Euler Eqs equilibrium points
% evaluate the stability switching the spinning axis direction and playing
% with different values of the off-axis initial angular velocities.

%% set initial conditions (switch i = x, y, z)
ix = 0.025; % kg m^2 minor axis (stable with no energy loss)
iy = 0.055; % kg m^2 intermediate axis (never stable)
iz = 0.070; % kg m^2 major axis (always stable)

wi0 = 2*pi; % rad s^-1, i = spinning axis

for i = 1:1:3
    w0 = [0.01 0.01 0.01];
    w0(i) = wi0;

    wx0 = w0(1);
    wy0 = w0(2);
    wz0 = w0(3);

    %% for verification
    wxdot0 = ((iy - iz) / ix) * wy0 * wz0;
    wydot0 = ((iz - ix) / iy) * wz0 * wx0;
    wzdot0 = ((ix - iy) / iz) * wx0 * wy0;
    
    %% get simulation outputs
    simout = sim('lab_3_task_3_simulink');
    wx = simout.wx;
    wy = simout.wy;
    wz = simout.wz;
    wxdot = simout.wxdot;
    wydot = simout.wydot;
    wzdot = simout.wzdot;
    
    %% plot data
    figure()
    plot(wxdot, "g")
    xlabel("Time (s)")
    ylabel("Angular Acceleration (rad s^-2)")
    title("Rotational Motion for a 3U Cubesat")
    grid on
    hold on
    plot(wydot, "b")
    plot(wzdot, "r")
    legend("wxdot", "wydot", "wzdot")
    hold off
    
    %% verify initial conditions
    disp("=== CONTROL INTIAL CONDITIONS ===")
    err = 1 - (wxdot.Data(1) / wxdot0);
    fprintf("wxdot error: %d CALC %d SIM: %d\n", err, wxdot0, wxdot.Data(1))
    err = 1 - (wydot.Data(1) / wydot0);
    fprintf("wydot error: %d CALC %d SIM: %d\n", err, wydot0, wydot.Data(1))
    err = 1 - (wzdot.Data(1) / wzdot0);
    fprintf("wzdot error: %d CALC %d SIM: %d\n", err, wzdot0, wzdot.Data(1))
    disp("=== CONTROL INTIAL CONDITIONS ===")
end