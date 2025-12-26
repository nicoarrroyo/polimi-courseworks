function B_N = get_mag_order_5_vectorised(R_planet, r_sat, latitude, longitude, G, H, order_of_model)
    % Inputs:
    %   latitude, longitude: vectors [deg]
    %   r_sat: vector [km]
    %   others: scalars
    
    latitude = latitude(:)';
    longitude = longitude(:)';
    r_sat = r_sat(:)'; % force row vectors
    N = length(latitude);

    colatitude = deg2rad(90 - latitude); 
    longitude_rad = deg2rad(longitude);

    Br = zeros(1, N); 
    Bt = zeros(1, N); 
    Bp = zeros(1, N);

    % Pre-calculate derivative step for P_nm (VECTORISED)
    d_theta = 1e-5;
    cos_theta = cos(colatitude);
    cos_theta_plus = cos(colatitude + d_theta);
    sin_colat = sin(colatitude);

    for n = 1:order_of_model
        % Distance ratio
        ratio = (R_planet ./ r_sat).^(n+2);
        
        % VECTORISED legendre
        P_n = legendre(n, cos_theta, "sch"); 
        P_n_plus = legendre(n, cos_theta_plus, "sch");
        dP_n = (P_n_plus - P_n) / d_theta;

        for m = 0:n
            g_nm = G(n, m+1);
            h_nm = H(n, m+1);
            
            % Get rows corresponding to order m (1 x N vectors)
            P_nm = P_n(m+1, :);
            dP_nm = dP_n(m+1, :);
            
            % VECTORISED Trig terms
            cos_m_phi = cos(m * longitude_rad);
            sin_m_phi = sin(m * longitude_rad);
            
            % Common term
            gh_term = (g_nm * cos_m_phi + h_nm * sin_m_phi);
            
            % calculate each component
            Br = Br + (ratio .* (n + 1) .* gh_term .* P_nm);
            Bt = Bt - (ratio .* gh_term .* dP_nm);
            Bp = Bp + (ratio .* m .* (-g_nm * sin_m_phi + h_nm * cos_m_phi) .* P_nm);
        end
    end

    % Apply the 1/sin(colatitude) factor to B_phi
    Bp = Bp .* (-1 ./ sin_colat);

    % --- LVLH -> ECI Transformation (Vectorized) ---
    % Construct rotation elements manually to allow vectorization
    % A_NLVLH = [ row1; row2; row3 ]
    
    sin_col = sin_colat; 
    cos_col = cos_theta;
    sin_lon = sin(longitude_rad);
    cos_lon = cos(longitude_rad);

    % Performing Matrix * Vector for N items:
    % B_N = A * B_LVLH
    
    % Row 1 of A * B_LVLH
    Bx = (sin_col .* cos_lon) .* Br + ...
         (cos_col .* cos_lon) .* Bt + ...
         (-sin_lon)           .* Bp;
         
    % Row 2 of A * B_LVLH
    By = (sin_col .* sin_lon) .* Br + ...
         (cos_col .* sin_lon) .* Bt + ...
         (cos_lon)            .* Bp;
         
    % Row 3 of A * B_LVLH
    Bz = (cos_col)            .* Br + ...
         (-sin_col)           .* Bt + ...
         0; % 0 * Bp

    % Output as N x 3 matrix
    B_N = [Bx', By', Bz'];
end