function acc_rsw = acc_pert_j2(t, s, params)
    % extract position and velocity vectors + get cartesian coordinates
    r_vec = s(1:3, :); v_vec = s(4:6, :);
    x = r_vec(1); y = r_vec(2); z = r_vec(3);
    r_norm = norm(r_vec);
    
    % 1. Compute J2 acceleration in Cartesian (ECI) 
    factor = (3/2) * params.J2 * params.mu * params.r_planet^2 / r_norm^4;
    aj2_x = factor * (x / r_norm) * (5 * (z^2 / r_norm^2) - 1);
    aj2_y = factor * (y / r_norm) * (5 * (z^2 / r_norm^2) - 1);
    aj2_z = factor * (z / r_norm) * (5 * (z^2 / r_norm^2) - 3);
    
    acc_eci = [aj2_x; aj2_y; aj2_z];
    
    % 2. Rotate to RSW Frame
    r_hat = r_vec / r_norm;
    
    h_vec = cross(r_vec, v_vec);
    h_norm = norm(h_vec);
    w_hat = h_vec / h_norm;
    
    s_hat = cross(w_hat, r_hat);
    
    % Rotation Matrix (ECI -> RSW)
    R_eci2rsw = [r_hat'; s_hat'; w_hat'];
    
    % Output acceleration in RSW
    acc_rsw = R_eci2rsw * acc_eci;
end
