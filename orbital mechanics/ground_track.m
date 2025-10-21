function [ long, lat ] = ground_track( y, t, w_E )
%% INPUTS
    % y, state of orbit at initial time (r, v)
    % theta_g0 [rad], longitude at greenwich meridian at initial time (for
    % simplicity this was set to 0)
    % t, vector of times at which the ground track will be computed
    % other inputs that you consider useful (rotation of earth w_E, earth
    % gravitational parameter, etc. (initial time set to 0))

%% OUTPUTS
    % alpha [rad], right ascensions in ECIEq
    % delta [rad], declination in ECIEq
    % long [rad], longitude phi = delta
    % lat [rad], latitude lambda = alpha - theta_G

%% state vector
r = y(1);

%% declination
delta = asin(r(3) / norm(r));

%% right ascension
if r(2) / norm(r) > 0
    alpha = acos(r(1) / (norm(r) * cos(delta)));
elseif r(2) / norm(r) <= 0
    alpha = 2*pi - acos(r(1) / (norm(r) * cos(delta)));
end

%% current longitude at greenwich meridian
theta_g = w_E * (t);

%% longitude
long = alpha - theta_g;

%% latitude
lat = delta;

