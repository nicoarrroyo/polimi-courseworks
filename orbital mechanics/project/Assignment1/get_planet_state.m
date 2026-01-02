function [r, v] = get_planet_state(mjd2000, body_id)
% GET_PLANET_STATE Get position and velocity of a planet
%
% Inputs:
%   mjd2000  - Modified Julian Date 2000 [days]
%   body_id  - Celestial body identifier
%
% Outputs:
%   r - Position vector [km] (column vector)
%   v - Velocity vector [km s^-1] (column vector)

    [kep, ~] = uplanet(mjd2000, body_id);

    mu = astroConstants(4);

    [r, v] = kep2car(...
        kep(1), ...  % Semi-major axis
        kep(2), ...  % Eccentricity
        kep(3), ...  % Inclination
        kep(4), ...  % RAAN
        kep(5), ...  % Argument of periapsis
        kep(6), ...  % True anomaly
        mu, ...
        "radians");
    
    r = r'; v = v'; % transpose for Lambert function
end
