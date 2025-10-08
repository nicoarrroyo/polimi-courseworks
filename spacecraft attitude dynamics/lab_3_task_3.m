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

I = diag([ix iy iz]);

wi0 = 2*pi; % rad s^-1, i = spinning axis

for i = 1:1:3
    w0 = [0.01 0.01 0.01];
    w0(i) = wi0;
    
    %% get simulation outputs
    simout = sim('lab_3_task_3_simulink');
    w = simout.w.Data;
    wxsim = w(:, 1);
    wysim = w(:, 2);
    wzsim = w(:, 3);

    wdot = simout.wdot.Data;
    wxdot = wdot(:, 1);
    wydot = wdot(:, 2);
    wzdot = wdot(:, 3);
    
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
end
