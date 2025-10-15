clear; close all; clc;

% Exercise 1 - 2-Body Problem
% Numerically integrate a Keplerian orbit (two-body problem)
%% %%% 1. Implement the code to propogate the orbit %%%%%%%%%%%%%%%%%%%%%%%
% a. identify states of the system and physical parameters involved
% b. write second-order ODE describing dynamics
% c. reduce problem to first-order ODE
% d. implement odefun function
% e. write a main script to numerically integrate the system
%     - choose one of MATLAB's solvers

mu_E = astroConstants(13); % [km^3 s^-2]

% initial conditions
r0 = [26578.137 0 0]; % [km]
v0 = [0 2.221 3.137]; % [km s^-1]
y0 = [r0 v0];

% time span
a = 1/(2/norm(r0) - dot(v0,v0)/mu_E); % [km]
T = 2*pi*sqrt(a^3/mu_E); % [s]
tspan = linspace(0, 5*T, 5000);

% set ODE solver conditions
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% integrate
[T, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y0, options);
r = Y(:, 1:3);
v = Y(:, 4:6);

%% %%% 2. Propagate the orbit for different initial  %%%%%%%%%%%%%%%%%%%%%%
% %%% conditions and primary attracting body %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a. for a closed orbit, initial conditions should correspond to a 
% negative specific energy. 
% did this; removed to save space

%% %%% 3. Analyse the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate initial specific energy
specific_energy_0 = 0.5 * norm(v0)^2 - mu_E / norm(r0);

% pre-allocate array sizes for efficiency
h = zeros(size(r));
e = zeros(size(r));
v_radial = zeros(length(T), 1);
v_transverse = zeros(length(T), 1);
e_dot_h = zeros(length(T), 1);
specific_energy = zeros(length(T), 1);

% calculate all the necessary variables
for i = 1:length(T)
    h(i, :) = cross(r(i, :), v(i, :));
    e(i, :) = (1 / mu_E) * cross(v(i, :), h(i, :)) - r(i, :) / norm(r(i, :));

    u_radial = r(i, :) / norm(r(i, :)); % radial unit vector
    u_oop = h(i, :) / norm(h(i, :)); % out-of-plane unit vector
    u_transverse = cross(u_oop, u_radial); % transverse unit vector
    v_radial(i) = dot(v(i, :), u_radial); % radial velocity
    v_transverse(i) = dot(v(i, :), u_transverse); % transverse velocity

    e_dot_h(i) = dot(e(i, :), h(i, :));

    specific_energy(i) = 0.5 * norm(v(i, :))^2 - mu_E / norm(r(i, :));
end

% %%% a. plot the orbit over 1 period T %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get texture of earth
earth_img = imread("EarthTexture.jpg");
[x_earth, y_earth, z_earth] = sphere(50);
radius = 6731;
x_earth = radius * x_earth;
y_earth = radius * y_earth;
z_earth = -radius * z_earth;

% plot
figure("Name", "Orbit Plot");
plot3(r(:, 1), r(:, 2), r(:, 3));
hold on
surface(x_earth, y_earth, z_earth, "FaceColor", "texturemap", ...
    "CData", earth_img, "EdgeColor", "none")
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Two-body problem orbit");
axis equal; grid on;
hold off

% %%% b. plot the components, norm of h, e, over 5 periods. they should be 
% %%% constant in magnitude and direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% angular momentum
figure("Name", "Angular Momentum")



% %%% c. check that h and e remain perpendicular by plotting the error %%%%


% %%% d. plot specific energy over 5 periods. should be constant in time %%


% %%% e. plot the evolution of vr and vθ (vt) over 5 periods %%%%%%%%%%%%%%




