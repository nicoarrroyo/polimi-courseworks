function B_LVLH = get_mag_order_5( R_planet, r_sat, latitude, longitude )

colatitude = deg2rad(latitude-90); longitude = deg2rad(longitude);

order_of_model = 5;
G = zeros(order_of_model, order_of_model + 1);
H = zeros(order_of_model, order_of_model + 1);

G(5, 1) = -234.42;
G(5, 2) = 47.52;    H(5, 2) = 47.52;
G(5, 3) = 208.36;   H(5, 3) = 208.36;
G(5, 4) = -121.43;  H(5, 4) = -121.43;
G(5, 5) = 32.09;    H(5, 5) = 32.09;
G(5, 6) = 13.98;    H(5, 6) = 99.14;

% Pre-calculate derivative step for P_nm
d_theta = 1e-5;
theta_plus = colatitude + d_theta;

Br = 0;
Bt = 0; % B_theta
Bp = 0; % B_phi

for n = 1:order_of_model
    % Calculate distance ratio (Re/r)^(n+2)
    ratio = (R_planet / r_sat)^(n+2);
    
    % Get Schmidt Semi-Normalized Legendre Polynomials for this n
    % P_n returns a vector of P_n^m for m = 0 to n
    P_n = legendre(n, cos(colatitude), "sch"); 
    
    % Calculate derivative using finite difference (robust method)
    P_n_plus = legendre(n, cos(theta_plus), "sch");
    dP_n = (P_n_plus - P_n) / d_theta;

    for m = 0:n
        % Extract coefficients (MATLAB indices start at 1)
        g_nm = G(n, m+1);
        h_nm = H(n, m+1);
        
        % Get specific Legendre value and derivative for this m
        P_nm = P_n(m+1);
        dP_nm = dP_n(m+1);
        
        % Common trig terms
        cos_m_phi = cos(m * longitude);
        sin_m_phi = sin(m * longitude);
        
        % --- Radial Component (Br) ---
        % Formula: Sum [ (n+1) * (g*cos + h*sin) * P_nm ]
        term_r = (n + 1) * (g_nm * cos_m_phi + h_nm * sin_m_phi) * P_nm;
        Br = Br + (ratio * term_r);
        
        % --- Co-elevation/colatitude Component (B_colatitude) ---
        % Formula: - Sum [ (g*cos + h*sin) * dP_nm/dcolatitude ]
        term_t = (g_nm * cos_m_phi + h_nm * sin_m_phi) * dP_nm;
        Bt = Bt - (ratio * term_t);
        
        % --- Azimuth/Phi Component (B_phi) ---
        % Formula: (-1/sin(colatitude)) * Sum [ m * (-g*sin + h*cos) * P_nm ]
        term_p = m * (-g_nm * sin_m_phi + h_nm * cos_m_phi) * P_nm;
        Bp = Bp + (ratio * term_p);
    end
end

% Apply the 1/sin(colatitude) factor to B_phi outside the summation
Bp = Bp * (-1 / sin(colatitude));

B_LVLH = [Br; Bt; Bp;];
end