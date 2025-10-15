clear; close all; clc;

%% Exercise 2 - Perturbed 2-Body Problem
% Earth orbit propogation with J2
J2 = astroConstants(9); % second zonal harmonic
mu_E = astroConstants(13); % [km^3 s^-2]
R_e = astroConstants(23); % earth radius [km]

% initial conditions
r0 = [6495 -970 -3622]; % [km]
v0 = [4.752 2.130 7.950]; % [km s^-1]
y0 = [r0 v0];

% time span
a = 1/(2/norm(r0) - dot(v0,v0)/mu_E); % [km]
tspan = linspace(0, 365*24*3600, 5000);

% set ODE solver conditions
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% integrate for unperturbed 2bp
[T_2bp, Y_2bp] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y0, options);
r_2bp = Y_2bp(:, 1:3);
v_2bp = Y_2bp(:, 4:6);

% integrate for j2-perturbed 2bp
[T_j2, Y_j2] = ode113(@(t,y) ode_2bp_j2(t,y,mu_E,J2,R_e), tspan, y0, options);
r_j2 = Y_j2(:, 1:3);
v_j2 = Y_j2(:, 4:6);

%% %%% 3. Analyse the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate initial specific energy
specific_energy_0 = 0.5 * norm(v0)^2 - mu_E / norm(r0);

% pre-allocate array sizes for efficiency
h = zeros(size(r_j2));
e = zeros(size(r_j2));
v_radial = zeros(length(T_j2), 1);
v_transverse = zeros(length(T_j2), 1);
e_dot_h = zeros(length(T_j2), 1);
specific_energy = zeros(length(T_j2), 1);

% calculate all the necessary variables
for i = 1:length(T_j2)
    h(i, :) = cross(r_j2(i, :), v_j2(i, :));
    e(i, :) = (1 / mu_E) * cross(v_j2(i, :), h(i, :)) - r_j2(i, :) / norm(r_j2(i, :));

    u_radial = r_j2(i, :) / norm(r_j2(i, :)); % radial unit vector
    u_oop = h(i, :) / norm(h(i, :)); % out-of-plane unit vector
    u_transverse = cross(u_oop, u_radial); % transverse unit vector
    v_radial(i) = dot(v_j2(i, :), u_radial); % radial velocity
    v_transverse(i) = dot(v_j2(i, :), u_transverse); % transverse velocity

    e_dot_h(i) = dot(e(i, :), h(i, :));

    specific_energy(i) = 0.5 * norm(v_j2(i, :))^2 - mu_E / norm(r_j2(i, :));
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
plot3(r_j2(:, 1), r_j2(:, 2), r_j2(:, 3));
hold on
patch(r_j2(:, 1), r_j2(:, 2), r_j2(:, 3), T_j2, "FaceColor", "none", "EdgeColor", "interp")
colorbar;
surface(x_earth, y_earth, z_earth, "FaceColor", "texturemap", ...
    "CData", earth_img, "EdgeColor", "none")
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
title("Two-body problem orbit");
axis equal; grid on;
hold off

% %%% b. plot the components, norm of h, e, over 5 periods. they should be 
% %%% constant in magnitude and direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% angular momentum
figure("Name", "Specific Angular Momentum")
plot(T_j2, h, "--");
hold on
plot(T_j2, vecnorm(h, 2, 2), "--");
xlabel("Time (s)"); ylabel("Specific Angular Momentum (km^2 s^-^1)");
title("Specific Angular Momentum Components and Norm");
legend("hx", "hy", "hz", "||h||");
ylim(ylim*1.1);
hold off

% eccentricity
figure("Name", "Eccentricity")
plot(T_j2, e, "--");
hold on
plot(T_j2, vecnorm(e, 2, 2), "--");
xlabel("Time (s)"); ylabel("Eccentricity (-)");
title("Eccentricity Components and Norm");
legend("ex", "ey", "ez", "||e||");
ylim(ylim*1.1);
hold off

% %%% c. check that h and e remain perpendicular by plotting the error %%%%
figure("Name", "Orthogonality Error of e and h")
plot(T_j2, e_dot_h);
xlabel("Time (s)"); ylabel("e dot h (km^2 s^-^1)");
title("Dot Product of e and h");

% %%% d. plot specific energy over 5 periods. should be constant in time %%
figure("Name", "Specific Energy")
plot(T_j2, specific_energy, "LineWidth", 0.8);
hold on
yline(specific_energy_0, "--")
xlabel("Time (s)"); ylabel("Specific Energy (km^2 s^-^2)");
title("Specific Energy Over 5 Periods");
legend("Calculated Specific Energy", "Initial Specific Energy");
hold off

% %%% e. plot the evolution of vr and vθ (vt) over 5 periods %%%%%%%%%%%%%%
figure("Name", "Radial and Transverse Velocity")
plot(T_j2, v_radial);
hold on
plot(T_j2, v_transverse);
xlabel("Time (s)"); ylabel("Velocity (km s^-^1)");
title("Radial and Transverse Velocity Over 5 Periods");
ylim(1.1*ylim); legend("Radial Velocity", "Transverse Velocity");
hold off
