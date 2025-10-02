%% Exercise 1 - 2-Body Problem
% Numerically integrate a Keplerian orbit (two-body problem)
% 1. identify states of the system and physical parameters involved
% 2. write second-order ODE describing dynamics
% 3. reduce problem to first-order ODE
% 4. implement odefun function
% 5. write a main script to numerically integrate the system
%     - choose one of MATLAB's solvers

% tspan = [tstart tend];
% y0 = y(tstart);
% [t, y] = odeXX(odefun, tspan, y0, options);



% initial conditions
r0 = [26578.137 0 0]; % [km]
v0 = [0 2.221 3.137]; % [km s^-1]
y0 = [r0 v0];

% time span


