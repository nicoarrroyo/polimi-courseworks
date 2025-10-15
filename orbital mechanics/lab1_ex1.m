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
tspan = linspace(0, 2*T, 1000);

% set ODE solver conditions
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% integrate
[T, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y0, options);

% get texture of earth
earth_img = imread("EarthTexture.jpg");
[x_earth, y_earth, z_earth] = sphere(50);
radius = 6731;
x_earth = radius * x_earth;
y_earth = radius * y_earth;
z_earth = -radius * z_earth;

% plot
figure()
plot3(Y(:,1), Y(:,2), Y(:,3), "-")
hold on
surface(x_earth, y_earth, z_earth, "FaceColor", "texturemap", ...
    "CData", earth_img, "EdgeColor", "none")
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Two-body problem orbit");
axis equal;
grid on;
hold off

%% %%% 2. Propagate the orbit for different initial  %%%%%%%%%%%%%%%%%%%%%%
% %%% conditions and primary attracting body %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a. for a closed orbit, initial conditions should correspond to a 
% negative specific energy. 
% did this; removed to save space

%% %%% 3. Analyse the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a. plot the orbit over 1 period T
figure()
plot3()

% b. plot the components, norm of h, e, over 5 periods. they should be 
% constant in magnitude and direction. 


% c. check that h and e remain perpendicular by plotting the error


% d. plot specific energy over 5 periods. should be constant in time. 


% e. plot the evolution of vr and vθ (vt) over 5 periods. 

% Calculate specific energy
specificEnergy = 0.5 * norm(v0)^2 - mu_E / norm(r0);




