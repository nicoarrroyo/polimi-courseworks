function [r0, v0] = kep2car( a, e, i, Omega, omega, TA)
    % ----- INPUTS -----
    % a: semi-major axis
    % e: eccentricity vector
    % i: inclination
    % RAAN: right ascension of ascending node
    % w (small omega): argument of perigee
    % TA: true anomaly
    % ------------------
    
    % ----- OUTPUTS -----
    % r0: cartesian position vector
    % v0: cartesian velocity vector
    % -------------------
    
    %% initial constants
    p = a * (1 - e^2);
    r = p / (1 + e * cosd(TA));
    h = cross(r, v);
    
    %% keplerian to perifocal (eph)
    % radius
    r_eph = [r * cosd(TA); ...
             r * sind(TA); ...
             0;]; % perifocal radius vector
    
    % velocity
    vr = (mu / h) * e * sind(TA); % radial velocity component
    vt = (mu / h) * (1 + e * cosd(TA)); % transverse velocity component
    v_eph = [vr * cosd(TA) - vt * sind(TA); ...
             vr * sind(TA) + vt * cosd(TA); ...
             0;]; % perifocal velocity vector
    
    %% perifocal to cartesian
    % use a 3-1-3 rotation matrix. in order:
        % 1. rotate around Z by the argument of perigee (omega)
        % 2. rotate around X by the inclination (i)
        % 3. rotate around Z by the RAAN (Omega)
    R3_omega = [cosd(omega) -sind(omega) 0;
                sind(omega) cosd(omega) 0;
                0 0 1];
    
    R1_i = [1 0 0;
            0 cosd(i) -sind(i);
            0 sind(i) cosd(i)];
    
    R3_Omega = [cosd(Omega) -sind(Omega) 0;
                sind(Omega) cosd(Omega) 0;
                0 0 1];
    
    r0 = R3_omega * R1_i * R3_Omega * r_eph;
    v0 = R3_omega * R1_i * R3_Omega * v_eph;
end
