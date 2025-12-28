function [dv1, dv2] = constrained_transfer_cost(...
    mjd2000_departure, ...
    mjd2000_arrival, ...
    body_id1, ...
    body_id2, ...
    dv_launch_constraint)
% TRANSFER_COST Calculate total delta-v for interplanetary transfer
%
%   dv_total = TRANSFER_COST(mjd2000_departure, mjd2000_arrival, body_id1, body_id2)
%   computes the total delta-v required for a Lambert transfer between two
%   celestial bodies.
%
% Inputs:
%   mjd2000_departure - Departure time [days from J2000]
%   mjd2000_arrival   - Arrival time [days from J2000]
%   body_id1          - Departure body ID (see below)
%   body_id2          - Arrival body ID
%
% Output:
%   dv_total - Total delta-v required [km s^-1]
%
% Body IDs:
%   1: Mercury,  2: Venus,    3: Earth,   4: Mars,    5: Jupiter
%   6: Saturn,   7: Uranus,   8: Neptune, 9: Pluto,  10: Sun

    % Constants
    mu_sun = astroConstants(4); % Sun's gravitational parameter [km^3 s^-2]

    % Convert time to seconds
    t1 = mjd2000_departure * 24 * 3600;
    t2 = mjd2000_arrival * 24 * 3600;
    tof = t2 - t1;

    % Get planetary states at departure and arrival
    [r1, v1] = get_planet_state(mjd2000_departure, body_id1, mu_sun);
    [r2, v2] = get_planet_state(mjd2000_arrival, body_id2, mu_sun);

    % Solve Lambert problem
    [~, ~, ~, ERROR, v_transfer1, v_transfer2, ~, ~] = ...
        lambertMR(r1, r2, tof, mu_sun, 0, 0, 0, 0);

    % Calculate delta-v maneuvers
    dv1 = abs(norm(v_transfer1 - v1));
    if ERROR == 0 && abs(dv1) < abs(dv_launch_constraint)
        dv2 = abs(norm(v_transfer2 - v2));
    else
       dv1 = NaN; dv2 = NaN;
    end
end

