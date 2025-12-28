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

    % backward compatibility safeguard to vectorisation update
    if size(r, 2) ~= 3 && size(r, 1) == 3
        r = r';
    end
    if size(v, 2) ~= 3 && size(v, 1) == 3
        v = v';
    end

    % 1. Calculate magnitudes
    r_mag = sqrt(sum(r.^2, 2));
    v_mag = sqrt(sum(v.^2, 2));

    % 2. Specific Angular Momentum (h)
    % Eq 1.12: h = r x v 
    h_vec = cross(r, v, 2);
    h = sqrt(sum(h_vec.^2, 2));

    % 3. Node Vector (n)
    % The line of nodes is the cross product of the Z-axis (K) and h
    % manually construct n_vec
    %K = [0, 0, 1];
    n_vec = [-h_vec(:, 2), h_vec(:, 1), zeros(size(h_vec, 1), 1)];
    n = sqrt(sum(n_vec.^2, 2));

    % 4. Eccentricity Vector (e)
    % Eq 1.21: e = (v x h)/mu - r/r 
    e_vec = (cross(v, h_vec, 2) ./ mu) - (r ./ r_mag);
    e = sqrt(sum(e_vec.^2, 2));

    % 5. Specific Mechanical Energy (E) and Semi-major Axis (a)
    % Eq 1.57: E = v^2/2 - mu/r = -mu/2a
    E = (v_mag.^2 ./ 2) - (mu ./ r_mag);
    
    if abs(E) > 1e-10 % Check for non-parabolic orbit
        a = -mu ./ (2 .* E);
    else
        a = inf; % Parabolic trajectory
    end

    % 6. Inclination (i)
    % Angle between K unit vector and h vector
    i = acos(h_vec(:, 3) ./ h);

    % 7. Right Ascension of the Ascending Node (Omega)
    % Angle between I unit vector ([1,0,0]) and n vector
    Omega = zeros(size(r,1), 1);

    non_eq = n > 1e-10;

    if any(non_eq)
        Omega(non_eq) = acos(n_vec(non_eq, 1) ./ n(non_eq));

        idx_quad = non_eq & (n_vec(:, 2) < 0);
        Omega(idx_quad) = 2*pi - Omega(idx_quad);
    end

    % 8. Argument of Periapsis (omega)
    omega = zeros(size(r, 1), 1);
    
    % Mask for defined omega (non-circular AND non-equatorial)
    valid_w = non_eq & (e > 1e-10);
    
    if any(valid_w)
        dot_ne = dot(n_vec(valid_w, :), e_vec(valid_w, :), 2);
        omega(valid_w) = acos(dot_ne ./ (n(valid_w) .* e(valid_w)));
        
        % Quadrant check: if ez < 0, omega > 180 deg
        idx_quad_w = valid_w & (e_vec(:, 3) < 0);
        omega(idx_quad_w) = 2*pi - omega(idx_quad_w);
    end

    % 9. True Anomaly (TA)
    TA = zeros(size(r, 1), 1);
    
    % Mask for non-circular orbits
    non_circ = e > 1e-10;
    
    if any(non_circ)
        dot_er = dot(e_vec(non_circ, :), r(non_circ, :), 2);
        TA(non_circ) = acos(dot_er ./ (e(non_circ) .* r_mag(non_circ)));
        
        % Quadrant check: if r.v < 0, TA > 180 deg
        dot_rv = dot(r(non_circ, :), v(non_circ, :), 2);
        idx_quad_ta = false(size(TA)); % Initialize full size mask
        
        % Map back to full array indices
        temp_indices = find(non_circ);
        sub_idx = dot_rv < 0;
        idx_quad_ta(temp_indices(sub_idx)) = true;
        
        TA(idx_quad_ta) = 2*pi - TA(idx_quad_ta);
    end
end
