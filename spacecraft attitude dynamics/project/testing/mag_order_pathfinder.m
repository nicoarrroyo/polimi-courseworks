clear; close all; clc; 

%% ========================================================================
%  ===== 1. ORGANIZE COEFFICIENTS INTO MATRICES =====
%  ========================================================================
% We create matrices G and H where the row index is n and column is m+1.
% Size: (max_order) x (max_order + 1)
order_of_model = 7;
G = zeros(order_of_model, order_of_model + 1);
H = zeros(order_of_model, order_of_model + 1);

% --- Order 1 ---
G(1, 1) = -29403.41;
G(1, 2) = -1451.37; H(1, 2) = 4653.35;

% --- Order 2 ---
G(2, 1) = -2499.78;
G(2, 2) = 2981.96;  H(2, 2) = -2991.72;
G(2, 3) = 1676.85;  H(2, 3) = -734.62;

% --- Order 3 ---
G(3, 1) = 1363.00;
G(3, 2) = -2380.80; H(3, 2) = -81.96;
G(3, 3) = 1236.06;  H(3, 3) = 241.80;
G(3, 4) = 525.60;   H(3, 4) = -542.52;

% --- Order 4 ---
G(4, 1) = 902.82;
G(4, 2) = 809.47;   H(4, 2) = 282.10;
G(4, 3) = 86.18;    H(4, 3) = -158.50;
G(4, 4) = -309.47;  H(4, 4) = 199.75;
G(4, 5) = 47.44;    H(4, 5) = -350.30;

% --- Order 5 ---
G(5, 1) = -234.42;
G(5, 2) = 47.52;    H(5, 2) = 47.52;
G(5, 3) = 208.36;   H(5, 3) = 208.36;
G(5, 4) = -121.43;  H(5, 4) = -121.43;
G(5, 5) = 32.09;    H(5, 5) = 32.09;
G(5, 6) = 13.98;    H(5, 6) = 99.14;

% --- Order 6 ---
G(6, 1) = 65.97;
G(6, 2) = 65.56;    H(6, 2) = -19.22;
G(6, 3) = 72.96;    H(6, 3) = 25.02;
G(6, 4) = -121.57;  H(6, 4) = 52.76;
G(6, 5) = -36.06;   H(6, 5) = -64.40;
G(6, 6) = 13.60;    H(6, 6) = 8.96;
G(6, 7) = -64.80;   H(6, 7) = 68.04;

% --- Order 7 ---
G(7, 1) = 80.54;
G(7, 2) = -76.63;   H(7, 2) = -51.50;
G(7, 3) = -8.23;    H(7, 3) = -16.85;
G(7, 4) = 56.45;    H(7, 4) = 2.36;
G(7, 5) = 15.80;    H(7, 5) = 23.56;
G(7, 6) = 6.30;     H(7, 6) = -2.19;
G(7, 7) = -7.21;    H(7, 7) = -27.19;
G(7, 8) = 9.77;     H(7, 8) = -1.90;

%% ========================================================================
%  ===== 2. DEFINE ORBIT PARAMETERS =====
%  ========================================================================
Re = 6371.2; % Earth Radius (km) - matches DGRF unit consistency
r_sat = 630 + Re; % Satellite radius (km)
theta = deg2rad(98.7-90); % Co-latitude (0 = North Pole, 180 = South)
phi_lon = deg2rad(15); % Longitude (East) (0 = Greenwich Meridian)

Br_list = []; Bt_list = []; Bp_list = []; tot_list = [];

%% ========================================================================
%  ===== 3. CALCULATION LOOP =====
%  ========================================================================
for i = 1:order_of_model
    max_order = i;

    % Initialize magnetic field components (in nT)
    Br = 0;
    Bt = 0; % B_theta
    Bp = 0; % B_phi
    
    % Pre-calculate derivative step for P_nm
    d_theta = 1e-5;
    theta_plus = theta + d_theta;
    
    for n = 1:max_order
        % Calculate distance ratio (Re/r)^(n+2)
        ratio = (Re / r_sat)^(n+2);
        
        % Get Schmidt Semi-Normalized Legendre Polynomials for this n
        % P_n returns a vector of P_n^m for m = 0 to n
        P_n = legendre(n, cos(theta), "sch"); 
        
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
            cos_m_phi = cos(m * phi_lon);
            sin_m_phi = sin(m * phi_lon);
            
            % --- Radial Component (Br) ---
            % Formula: Sum [ (n+1) * (g*cos + h*sin) * P_nm ]
            term_r = (n + 1) * (g_nm * cos_m_phi + h_nm * sin_m_phi) * P_nm;
            Br = Br + (ratio * term_r);
            
            % --- Co-elevation/Theta Component (B_theta) ---
            % Formula: - Sum [ (g*cos + h*sin) * dP_nm/dtheta ]
            term_t = (g_nm * cos_m_phi + h_nm * sin_m_phi) * dP_nm;
            Bt = Bt - (ratio * term_t); % Note the negative sign
            
            % --- Azimuth/Phi Component (B_phi) ---
            % Formula: (-1/sin(theta)) * Sum [ m * (-g*sin + h*cos) * P_nm ]
            term_p = m * (-g_nm * sin_m_phi + h_nm * cos_m_phi) * P_nm;
            Bp = Bp + (ratio * term_p);
        end
    end
    
    % Apply the 1/sin(theta) factor to B_phi outside the summation
    Bp = Bp * (-1 / sin(theta));
    
    %% ====================================================================
    %  ===== 4. RESULTS AND CONVERSION =====
    %  ====================================================================
    % The result [Br; Bt; Bp] is in the LOCAL HORIZONTAL 
    % frame (NEC/NED depending on definition)
    % B_local = [Br; Bt; Bp];

    Br_list = [Br_list; Br;];
    Bt_list = [Bt_list; Bt;];
    Bp_list = [Bp_list; Bp;];

    tot = norm([Br Bt Bp]);
    tot_list = [tot_list; tot;];

    if length(Br_list) <= 1
        diff = [0; 0; 0];
        diff_tot = 0;
    else
        diff_r = 100 * (Br_list(end) - Br_list(end-1)) / Br_list(end-1);
        diff_t = 100 * (Bt_list(end) - Bt_list(end-1)) / Bt_list(end-1);
        diff_p = 100 * (Bp_list(end) - Bp_list(end-1)) / Bp_list(end-1);
        diff = [diff_r; diff_t; diff_p];

        diff_tot = 100 * (tot_list(end) - tot_list(end-1)) / tot_list(end-1);
    end
    
    fprintf('\nOrder: %.1f\n', max_order);
    fprintf("===========\n")
    fprintf('B_radial: %.2f\n', Br);
    fprintf("Percentage diff. from pervious order: %.2f\n", diff(1))
    fprintf('B_theta:  %.2f\n', Bt);
    fprintf("Percentage diff. from pervious order: %.2f\n", diff(2))
    fprintf('B_phi:    %.2f\n', Bp);
    fprintf("Percentage diff. from pervious order: %.2f\n", diff(3))
    fprintf('Total Intensity: %.2f\n', tot);
    fprintf("Percentage diff. from pervious order: %.2f\n", diff_tot)
end

%% ========================================================================
%  ===== 5. PLOTS =====
%  ========================================================================
figure("Name", "Radial Component (B_r)")
plot(Br_list)

figure("Name", "Co-Elevation Component (B_theta)")
plot(Bt_list)


figure("Name", "Azimuth component (B_phi)")
plot(Bp_list)


