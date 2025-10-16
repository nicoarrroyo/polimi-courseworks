clear; close all; clc;

%% Exercise 2 - Perturbed 2-Body Problem
% Earth orbit propogation with J2
J2 = astroConstants(9); % second zonal harmonic
mu_E = astroConstants(13); % [km^3 s^-2]
R_e = astroConstants(23); % earth radius [km]

% initial conditions
r0 = [6495 -970 -3622]; % initial position vector [km]
v0 = [4.752 2.130 7.950]; % initial velocity vector [km s^-1]
y0 = [r0 v0]; % initial state vector state vector

% others
a = 1/(2/norm(r0) - dot(v0,v0)/mu_E); % semi major axis [km]
tspan = linspace(0, 365*24*3600, 50000); % time span

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

% pre-allocate array sizes for efficiency unperturbed 2bp
h_2bp = zeros(size(r_2bp));
e_2bp = zeros(size(r_2bp));
v_radial_2bp = zeros(length(T_2bp), 1);
v_transverse_2bp = zeros(length(T_2bp), 1);
e_dot_h_2bp = zeros(length(T_2bp), 1);
specific_energy_2bp = zeros(length(T_2bp), 1);

% calculate all the necessary variables unperturbed 2bp
for i = 1:length(T_j2)
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
h = zeros(size(r_j2));
e = zeros(size(r_j2));
v_radial = zeros(length(T_j2), 1);
v_transverse = zeros(length(T_j2), 1);
e_dot_h = zeros(length(T_j2), 1);
specific_energy = zeros(length(T_j2), 1);

% calculate all the necessary variables j2
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
plot3(r_2bp(:, 1), r_2bp(:, 2), r_2bp(:, 3), "r", "LineWidth", 4);
hold on
patch(r_j2(:, 1), r_j2(:, 2), r_j2(:, 3), T_j2/(24*3600), "FaceColor", ...
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
plot(T_2bp/(24*3600), h_2bp, "--");
hold on
plot(T_2bp/(24*3600), vecnorm(h_2bp, 2, 2), "--");
plot(T_j2/(24*3600), h, "--");
plot(T_j2/(24*3600), vecnorm(h, 2, 2), "--");
xlabel("Time (days)"); ylabel("Specific Angular Momentum (km^2 s^-^1)");
title("Specific Angular Momentum Components and Norm");
legend("hx 2BP", "hy 2BP", "hz 2BP", "||h|| 2BP", ...
    "hx J2", "hy J2", "hz J2", "||h|| J2");
ylim(ylim*1.1);
hold off

% eccentricity
figure("Name", "Eccentricity")
plot(T_2bp/(24*3600), e_2bp, "--");
hold on
plot(T_2bp/(24*3600), vecnorm(e_2bp, 2, 2), "--");
plot(T_j2/(24*3600), e, "o");
plot(T_j2/(24*3600), vecnorm(e, 2, 2), "--");
xlabel("Time (days)"); ylabel("Eccentricity (-)");
title("Eccentricity Components and Norm");
legend("ex 2BP", "ey 2BP", "ez 2BP", "||e|| 2BP", ...
    "ex J2", "ey J2", "ez J2", "||e|| J2");
ylim(ylim*1.1);
hold off

% %%% c. check that h and e remain perpendicular by plotting the error %%%%
figure("Name", "Orthogonality Error of e and h")
plot(T_2bp/(24*3600), e_dot_h_2bp);
hold on
plot(T_j2/(24*3600), e_dot_h);
xlabel("Time (days)"); ylabel("e dot h (km^2 s^-^1)");
legend("e dot h 2BP", "e dot h J2");
title("Dot Product of e and h");

% %%% d. plot specific energy over 5 periods. should be constant in time %%
figure("Name", "Specific Energy")
plot(T_2bp/(24*3600), specific_energy_2bp, "LineWidth", 0.8);
hold on
plot(T_j2/(24*3600), specific_energy, "LineWidth", 0.8);
yline(specific_energy_0, "--")
xlabel("Time (days)"); ylabel("Specific Energy (km^2 s^-^2)");
title("Specific Energy Over 5 Periods");
legend("Calculated Specific Energy 2BP", ...
    "Calculated Specific Energy J2", "Initial Specific Energy");
hold off

% %%% e. plot the evolution of vr and vθ (vt) over 5 periods %%%%%%%%%%%%%%
figure("Name", "Radial and Transverse Velocity")
plot(T_2bp/(24*3600), v_radial_2bp);
hold on
plot(T_2bp/(24*3600), v_transverse_2bp);
plot(T_j2/(24*3600), v_radial);
plot(T_j2/(24*3600), v_transverse);
xlabel("Time (days)"); ylabel("Velocity (km s^-^1)");
title("Radial and Transverse Velocity Over 5 Periods");
ylim(1.1*ylim); 
legend("Radial Velocity 2BP", "Transverse Velocity 2BP", ...
    "Radial Velocity J2", "Transverse Velocity J2");
hold off