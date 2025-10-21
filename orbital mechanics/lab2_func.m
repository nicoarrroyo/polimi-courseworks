function lab2_func( r0, v0, orbits, plot_title )
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