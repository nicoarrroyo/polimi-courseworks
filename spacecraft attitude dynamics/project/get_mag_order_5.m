function B_N = get_mag_order_5( R_planet, r_sat, latitude, longitude, G, H, order_of_model )

colatitude = deg2rad(90 - latitude); longitude = deg2rad(longitude);

% Pre-calculate derivative step for P_nm
d_theta = 1e-5;
theta_plus = colatitude + d_theta;
Br = 0; Bt = 0; Bp = 0; % B_radial, B_theta, B_phi

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
        % Extract coefficients
        g_nm = G(n, m+1);
        h_nm = H(n, m+1);
        
        % Get specific Legendre value and derivative for this m
        P_nm = P_n(m+1);
        dP_nm = dP_n(m+1);
        
        % Common trig terms
        cos_m_phi = cos(m * longitude);
        sin_m_phi = sin(m * longitude);
        
        % --- Radial Component (Br) ---
        Br = Br + (ratio * ((n + 1) * (g_nm * cos_m_phi + h_nm * sin_m_phi) * P_nm));
        
        % --- Co-elevation/colatitude Component (B_colatitude) ---
        Bt = Bt - (ratio * ((g_nm * cos_m_phi + h_nm * sin_m_phi) * dP_nm));
        
        % --- Azimuth/Phi Component (B_phi) ---
        Bp = Bp + (ratio * (m * (-g_nm * sin_m_phi + h_nm * cos_m_phi) * P_nm));
    end
end

% Apply the 1/sin(colatitude) factor to B_phi outside the summation
Bp = Bp * (-1 / sin(colatitude));
B_LVLH = [Br; Bt; Bp;];

% LVLH -> ECI
A_NLVLH = [ sin(colatitude)*cos(longitude),  cos(colatitude)*cos(longitude), -sin(longitude) ;
         sin(colatitude)*sin(longitude),  cos(colatitude)*sin(longitude),  cos(longitude) ;
         cos(colatitude),           -sin(colatitude),            0        ];

B_N = A_NLVLH * B_LVLH;

end
