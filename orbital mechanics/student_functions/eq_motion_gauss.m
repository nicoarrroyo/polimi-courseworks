function dkep_dt = eq_motion_gauss(t, kep, a_per, params)
    a = kep(1); e = kep(2); i = kep(3); Omega = kep(4); omega = kep(5); TA = kep(6);

    % Auxiliary parameters
    p = a * (1 - e^2); % semi-latus rectum
    h = sqrt(p * params.mu); % angular momentum
    r = p / (1 + e * cos(TA)); % position vector magntitude
    
    % 1. Convert current Keplerian state to Cartesian (r_vec, v_vec)
    [r_vec, v_vec] = kep2car( ...   % state vector [km, km s^-1]
        a, ...                      % a: semi major axis [km]
        e, ...                      % e: eccentricity [-]
        i, ...                      % i: inclination [rad]
        Omega, ...                  % Omega: right ascension of the ascending node [rad]
        omega, ...                  % omega: argument of perigee [rad]
        TA, ...                     % TA: true anomaly [rad]
        params.mu, ...
        "radians" ...               % must specify unit of angles if not deg
        );
    
    s_vec = [r_vec; v_vec;];

    % 2. Calculate Perturbing Acceleration in RSW frame
    % We pass the Cartesian state because J2 depends on Cartesian position (z component)
    % a_per is a generic acceleration perturbation function that can be
    % switched out with acc_pert_j2 or acc_pert_thirdbody, etc.
    acc_pert_vec = a_per(t, s_vec, params);
    ar = acc_pert_vec(1);
    as = acc_pert_vec(2);
    aw = acc_pert_vec(3);
    
    % 3. Gauss Equations (RSW Form) 
    da_dt = (2 * a^2 / h) * (e * sin(TA) * ar + (p / r) * as);

    de_dt = (1 / h) * (p * sin(TA) * ar + ( (p + r) * cos(TA) + r * e ) * as);

    di_dt = (r * cos(TA + omega) / h) * aw;

    dOmega_dt = (r * sin(TA + omega) / (h * sin(i))) * aw;

    domega_dt = (1 / (h * e)) * (-p * cos(TA) * ar + (p + r) * sin(TA) * as) ...
            - (r * sin(TA + omega) * cos(i) / (h * sin(i))) * aw;
            
    dTA_dt = h / r^2 + (1 / (e * h)) * (p * cos(TA) * ar - (p + r) * sin(TA) * as);

    dkep_dt = [da_dt; de_dt; di_dt; dOmega_dt; domega_dt; dTA_dt];
end