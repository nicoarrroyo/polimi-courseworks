function lab2_ground_track( r0, v0, orbits, plot_title )
%% configure paths
script_path = fileparts(mfilename("fullpath"));
backs = strfind(script_path, "\"); labs_dir = script_path(1:backs(end));
addpath([labs_dir '\student_functions']); addpath([labs_dir '\lib']);

%% constants / initial conditions
w_E = 15.04; % earth rotation velocity [deg hr^-1]
mu_E = astroConstants(13); % earth gravitational parameter [km^3 s^-2]
earth_img = imread("EarthTexture.jpg");

%% intial conditions
a = 1 / (2 / norm(r0) - dot(v0,v0) / mu_E); % semi-major axis [km]
y0 = [r0 v0];

%% time span
period = 2 * pi * sqrt(a^3 / mu_E); % orbital period [s]
tspan = linspace(0, orbits*period, 10000);

%% orbit propogation
% set solver conditions
options = odeset("RelTol", 1e-13, "AbsTol", 1e-14);

% integrate
[T, Y] = ode113(@(t,y) ode_2bp(t,y,mu_E), tspan, y0, options);

%% ground track calculation and plotting
% calculation
[long, lat] = ground_track(Y, T, w_E); % longitude and latitude [radians]
fprintf("FINAL LONGITUDE: %.4f rad, %.2f deg\n", long(end), rad2deg(long(end)));
fprintf("FINAL LATITUDE:  %.4f rad, %.2f deg\n", lat(end), rad2deg(lat(end)));

% plotting
figure("Name", plot_title)
image([-180 180], [90 -90], earth_img);

hold on
plot(rad2deg(long), rad2deg(lat), ".")
xline(rad2deg(long(1)), "r");
xline(rad2deg(long(end)), "g");

title(plot_title)
xlabel("longitude [deg]"); ylabel("latitude [deg]");
ylim([-90 90]); xlim([-180 180]);
legend("ground track", "initial long", "final long")
set(gca, "Ydir", "normal")
hold off