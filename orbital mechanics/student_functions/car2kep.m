function [a, e, i, Omega, omega, TA] = car2kep( r, v, mu )
    % ----- INPUTS -----
    % r0: cartesian position vector
    % v0: cartesian velocity vector
    % mu: gravitational parameter
    % ------------------

    % ----- OUTPUTS -----
    % a: semi-major axis
    % e: eccentricity scalar
    % i: inclination
    % RAAN: right ascension of ascending node
    % w (small omega): argument of perigee
    % TA: true anomaly
    % -------------------
    if nargin < 3
        mu = astroConstants(13); % earth gravitational parameter
    end

    % 1. Calculate magnitudes
    r_mag = norm(r);
    v_mag = norm(v);

    % 2. Specific Angular Momentum (h)
    % Eq 1.12: h = r x v 
    h_vec = cross(r, v);
    h = norm(h_vec);

    % 3. Node Vector (n)
    % The line of nodes is the cross product of the Z-axis (K) and h
    K = [0, 0, 1]; 
    n_vec = cross(K, h_vec);
    n = norm(n_vec);

    % 4. Eccentricity Vector (e)
    % Eq 1.21: e = (v x h)/mu - r/r 
    e_vec = (cross(v, h_vec) / mu) - (r / r_mag);
    e = norm(e_vec);

    % 5. Specific Mechanical Energy (E) and Semi-major Axis (a)
    % Eq 1.57: E = v^2/2 - mu/r = -mu/2a
    E = (v_mag^2 / 2) - (mu / r_mag);
    
    if abs(E) > 1e-10 % Check for non-parabolic orbit
        a = -mu / (2 * E);
    else
        a = inf; % Parabolic trajectory
    end

    % 6. Inclination (i)
    % Angle between K unit vector and h vector
    i = acos(h_vec(3) / h);

    % 7. Right Ascension of the Ascending Node (Omega)
    % Angle between I unit vector ([1,0,0]) and n vector
    if n ~= 0
        Omega = acos(n_vec(1) / n);
        % Quadrant check: if ny < 0, Omega > 180 deg
        if n_vec(2) < 0
            Omega = 2*pi - Omega;
        end
    else
        Omega = 0; % Undefined for equatorial orbits
    end

    % 8. Argument of Periapsis (omega)
    % Angle between node vector n and eccentricity vector e
    if n ~= 0 && e > 1e-10
        omega = acos(dot(n_vec, e_vec) / (n * e));
        % Quadrant check: if ez < 0, omega > 180 deg
        if e_vec(3) < 0
            omega = 2*pi - omega;
        end
    else
        omega = 0; % Undefined for circular or equatorial orbits
    end

    % 9. True Anomaly (theta)
    % Angle between eccentricity vector e and position vector r
    % Eq 1.26 definition context
    if e > 1e-10
        TA = acos(dot(e_vec, r) / (e * r_mag));
        % Quadrant check: if r.v < 0, the satellite is moving towards periapsis
        % (flying from apocentre to pericentre), so theta > 180 deg
        if dot(r, v) < 0
            TA = 2*pi - TA;
        end
    else
        TA = 0; % Undefined for circular orbits
    end

end