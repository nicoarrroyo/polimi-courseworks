function [ long, lat ] = ground_track( Y, T, w_E )
%% INPUTS
    % Y, state of orbit throughout timespan time (r, v)
    % theta_g0 [rad], longitude at greenwich meridian at initial time (for
    % simplicity this was set to 0)
    % T, vector of times at which the ground track will be computed
    % other inputs that you consider useful (rotation of earth w_E, earth
    % gravitational parameter, etc. (initial time set to 0))

%% OUTPUTS
    % alpha [rad], right ascensions in ECIEq
    % delta [rad], declination in ECIEq
    % long [rad], longitude phi = delta
    % lat [rad], latitude lambda = alpha - theta_G

%% state vector
r = Y(:, 1:3);
r_norm = vecnorm(r')';

%% declination
delta = asin(r(:, 3) ./ r_norm);

%% right ascension
% if r(:, 2) / r_norm > 0
%     alpha = acos(r(1) ./ (r_norm .* cos(delta)));
% elseif r(:, 2) / r_norm <= 0
%     alpha = 2*pi - acos(r(1) / (r_norm .* cos(delta)));
% else
%     disp("no value of alpha")
% end
alpha = atan2(r(:, 2), r(:, 1));

%% current axial longitude from greenwich meridian
theta_g = (deg2rad(w_E) / 3600) * T;

%% longitude
long = alpha - theta_g;
long = wrapToPi(long);

%% latitude
lat = delta;

