function dv_total = transfer_cost(mjd2000_departure, mjd2000_arrival, body_id1, body_id2)
% TRANSFER_COST Calculate total delta-v for interplanetary transfer
%
%   dv_total = TRANSFER_COST(mjd2000_departure, mjd2000_arrival, body_id1, body_id2)
%   computes the total delta-v required for a Lambert transfer between two
%   celestial bodies.
%
% Inputs:
%   mjd2000_departure - Departure time [days from J2000]
%   mjd2000_arrival   - Arrival time [days from J2000]
%   body_id1          - Departure body ID (1=Mercury, 2=Venus, 3=Earth, etc.)
%   body_id2          - Arrival body ID
%
% Output:
%   dv_total - Total delta-v required [km/s]
%
% Body IDs:
%   1: Mercury,  2: Venus,    3: Earth,   4: Mars,    5: Jupiter
%   6: Saturn,   7: Uranus,   8: Neptune, 9: Pluto,  10: Sun
    script_path = fileparts(mfilename('fullpath'));
    path_separators = strfind(script_path, '\');
    labs_dir = script_path(1:path_separators(end));

    addpath(fullfile(labs_dir, 'student_functions'));
    addpath(fullfile(labs_dir, 'lib'));
    addpath(fullfile(labs_dir, 'lib', 'timeConversion'));

    % Constants
    mu_sun = astroConstants(4); % Sun's gravitational parameter [km^3/s^2]
    SEC_PER_DAY = 86400;

    % Convert time to seconds
    t1 = mjd2000_departure * SEC_PER_DAY;
    t2 = mjd2000_arrival * SEC_PER_DAY;
    tof = t2 - t1;

    % Get planetary states at departure and arrival
    [r1, v1] = get_planet_state(mjd2000_departure, body_id1);
    [r2, v2] = get_planet_state(mjd2000_arrival, body_id2);

    % Solve Lambert problem
    [~, ~, ~, ~, v_transfer1, v_transfer2, ~, ~] = ...
        lambertMR(r1, r2, tof, mu_sun, 0, 0, 0, 0);

    % Calculate delta-v maneuvers
    dv1 = norm(v_transfer1 - v1);
    dv2 = norm(v_transfer2 - v2);
    dv_total = dv1 + dv2;
end

function [r, v] = get_planet_state(mjd2000, body_id)
% GET_PLANET_STATE Get position and velocity of a planet
%
% Inputs:
%   mjd2000  - Modified Julian Date 2000 [days]
%   body_id  - Celestial body identifier
%
% Outputs:
%   r - Position vector [km] (column vector)
%   v - Velocity vector [km/s] (column vector)

    % Get Keplerian elements
    [kep, ~] = uplanet(mjd2000, body_id);

    % Convert to Cartesian coordinates
    [r, v] = kep2car(...
        kep(1), ...  % Semi-major axis
        kep(2), ...  % Eccentricity
        kep(3), ...  % Inclination
        kep(4), ...  % RAAN
        kep(5), ...  % Argument of periapsis
        kep(6));     % True anomaly

    % Transpose to column vectors for Lambert function
    r = r';
    v = v';
end

