%% Exercise 1 - 2-Body Problem
% Numerically integrate a Keplerian orbit (two-body problem)
% 1. identify states of the system and physical parameters involved
% 2. write second-order ODE describing dynamics
% 3. reduce problem to first-order ODE
% 4. implement odefun function
% 5. write a main script to numerically integrate the system
%     - choose one of MATLAB's solvers

mu_E = astroConstants(13); % [km^3 s^-2]

% initial conditions
r0 = [26578.137 0 0]; % [km]
v0 = [0 2.221 3.137]; % [km s^-1]
y0 = [r0 v0];

% time span
a = 1/(2/norm(r0) - dot(v0,v0)/mu_E); % [km]
T = 2*pi*sqrt(a^3/mu_E); % [s]
tspan = linspace(0, 2*T, 1000);

% set ODE solver conditions
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% integrate
[T, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y0, options);

% plot
figure()
plot3(Y(:,1), Y(:,2), Y(:,3), "-")
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Two-body problem orbit");
axis equal;
grid on;

