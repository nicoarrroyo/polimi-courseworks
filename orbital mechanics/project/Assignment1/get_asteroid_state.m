function [r, v] = get_asteroid_state(mjd2000, body_id, mu)
% GET_PLANET_STATE Get position and velocity of a planet
%
% Inputs:
%   mjd2000  - Modified Julian Date 2000 [days]
%   body_id  - Asteroid identifier number
%   mu       - gravitational parameter of the body in question
%
% Outputs:
%   r - Position vector [km] (column vector)
%   v - Velocity vector [km s^-1] (column vector)

    [kep, ~] = uplanet(mjd2000, body_id);
    [kep, ~, M] = ephAsteroids(mjd2000, body_id);
    
    % --- convert from mean to eccentric anomaly anomaly ---
    % M = E - esinE -> solve iteratively
    tol = 10^-8;
    E = M; e = kep(2);
    ratio = 1;
    while abs(ratio) > tol
        f = E - e * sin(E) - M;
        f_prime = 1 - e * cos(E);

        ratio = f / f_prime;
        E = E - ratio;
    end

    % --- convert from eccentric anomaly to true anomaly ---
    % TA = 2arctan(tan(E/2) / sqrt((1-e)/(1+e)))
    I'M SO DUMB!! 
    it gives you a kep vector
    which has true anomaly! 
    no need to mess around with converting it from mean anomaly

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
