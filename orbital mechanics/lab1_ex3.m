clear; close all; clc;

%% Exercise 3 - Kepler's Equation
% Earth orbit propogation with J2
J2 = astroConstants(9); % second zonal harmonic
mu_E = astroConstants(13); % [km^3 s^-2]
R_e = astroConstants(23); % earth radius [km]

% initial conditions
r0 = [6495 -970 -3622]; % initial position vector [km]
v0 = [4.752 2.130 7.950]; % initial velocity vector [km s^-1]
y0 = [r0 v0]; % initial state vector state vector
E0 = 0; % initial eccentric anomaly [rad]
t0 = 0; % initial time [s]
e_values = [0 0.2 0.4 0.6 0.8 0.95];

% others
% a = 1/(2/norm(r0) - dot(v0,v0)/mu_E); % semi major axis [km]
a = 7000; % semi major axis [km]
k = 2; % periods of the orbit
N = 150; % points of orbit

tspan = linspace(0, 365*24*3600, N); % time span

% set ODE solver conditions
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% integrate for j2-perturbed 2bp
[T, Y] = ode113(@(t,y) ode_2bp_j2(t,y,mu_E,J2,R_e), tspan, y0, options);
r = Y(:, 1:3);
v = Y(:, 4:6);

%% %%% 3. Analyse the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate initial specific energy
specific_energy_0 = 0.5 * norm(v0)^2 - mu_E / norm(r0);

% pre-allocate array sizes for efficiency unperturbed 2bp
h_2bp = zeros(size(r_2bp));
e_2bp = zeros(size(r_2bp));
v_radial_2bp = zeros(length(T_2bp), 1);
v_transverse_2bp = zeros(length(T_2bp), 1);
e_dot_h_2bp = zeros(length(T_2bp), 1);
specific_energy_2bp = zeros(length(T_2bp), 1);

% calculate all the necessary variables unperturbed 2bp
for i = 1:length(T)
    h_2bp(i, :) = cross(r_2bp(i, :), v_2bp(i, :));
    e_2bp(i, :) = (1 / mu_E) * cross(v_2bp(i, :), h_2bp(i, :)) - ...
        r_2bp(i, :) / norm(r_2bp(i, :));

    u_radial = r_2bp(i, :) / norm(r_2bp(i, :)); % radial unit vector
    u_oop = h_2bp(i, :) / norm(h_2bp(i, :)); % out-of-plane unit vector
    u_transverse = cross(u_oop, u_radial); % transverse unit vector
    v_radial_2bp(i) = dot(v_2bp(i, :), u_radial); % radial velocity
    v_transverse_2bp(i) = dot(v_2bp(i, :), u_transverse); % transverse velocity

    e_dot_h_2bp(i) = dot(e_2bp(i, :), h_2bp(i, :));

    specific_energy_2bp(i) = 0.5 * norm(v_2bp(i, :))^2 - mu_E / norm(r_2bp(i, :));
end

% pre-allocate array sizes for efficiency j2
h = zeros(size(r));
e = zeros(size(r));
v_radial = zeros(length(T), 1);
v_transverse = zeros(length(T), 1);
e_dot_h = zeros(length(T), 1);
specific_energy = zeros(length(T), 1);

% calculate all the necessary variables j2
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
R_e = astroConstants(23); % earth radius [km]
x_earth = R_e * x_earth;
y_earth = R_e * y_earth;
z_earth = -R_e * z_earth;

% plot
figure("Name", "Orbit Plot");
plot3(r_2bp(:, 1), r_2bp(:, 2), r_2bp(:, 3), "r", "LineWidth", 4);
hold on
patch(r(:, 1), r(:, 2), r(:, 3), T/(24*3600), "FaceColor", ...
    "none", "EdgeColor", "interp")
colorbar;
surface(x_earth, y_earth, z_earth, "FaceColor", "texturemap", ...
    "CData", earth_img, "EdgeColor", "none")
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Two-body problem orbit");
axis equal; grid on;
legend("Unperturbed 2BP", "J2 Perturbed 2BP");
hold off

% %%% b. plot the components, norm of h, e, over 5 periods. they should be 
% %%% constant in magnitude and direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% angular momentum
figure("Name", "Specific Angular Momentum")
plot(T/(24*3600), h, "--");
hold on
plot(T/(24*3600), vecnorm(h, 2, 2), "--");
xlabel("Time (days)"); ylabel("Specific Angular Momentum (km^2 s^-^1)");
title("Specific Angular Momentum Components and Norm");
legend("hx", "hy", "hz", "||h||");
ylim(ylim*1.1);
hold off

% eccentricity
figure("Name", "Eccentricity")
plot(T/(24*3600), e, "o");
hold on
plot(T/(24*3600), vecnorm(e, 2, 2), "--");
xlabel("Time (days)"); ylabel("Eccentricity (-)");
title("Eccentricity Components and Norm");
legend("ex", "ey", "ez", "||e||");
ylim(ylim*1.1);
hold off

% %%% c. check that h and e remain perpendicular by plotting the error %%%%
figure("Name", "Orthogonality Error of e and h")
plot(T/(24*3600), e_dot_h);
xlabel("Time (days)"); ylabel("e dot h (km^2 s^-^1)");
legend("e dot h");
title("Dot Product of e and h");

% %%% d. plot specific energy over 5 periods. should be constant in time %%
figure("Name", "Specific Energy")
plot(T/(24*3600), specific_energy, "LineWidth", 0.8);
hold on
yline(specific_energy_0, "--")
xlabel("Time (days)"); ylabel("Specific Energy (km^2 s^-^2)");
title("Specific Energy Over 5 Periods");
legend("Calculated Specific Energy", "Initial Specific Energy");
hold off

% %%% e. plot the evolution of vr and vθ (vt) over 5 periods %%%%%%%%%%%%%%
figure("Name", "Radial and Transverse Velocity")
plot(T/(24*3600), v_radial);
hold on
plot(T/(24*3600), v_transverse);
xlabel("Time (days)"); ylabel("Velocity (km s^-^1)");
title("Radial and Transverse Velocity Over 5 Periods");
ylim(1.1*ylim); 
legend("Radial Velocity", "Transverse Velocity");
hold off