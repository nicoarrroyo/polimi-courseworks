clc, clearvars, close all

%% TEMPERATURE RANGE

% Needs the real operating values of all the onboard components
% Tmax is the minimum of the maximum temperatures
% Tmin is the maximum of the minumum temperatures

T_max_SC            = 273.15 + 50;
T_min_SC            = 273.15 - 20;

%% AREA AND ALPHA/EPS CALCULATIONS (TODO)

% A_tot calculation needed ---> turn into equivalent sphere
% Different absorbptivity and emissivity for the SC surface
% Weighted average needed

% Material properties
alpha_mli           = 0.39;
eps_mli             = 0.78;

% MIRAS arms (front: average LICEF + coating | back: coating)
A_arms              = 4.08;
alpha_MIRAS_front   = (0.30*0.1 + alpha_mli*0.9);
eps_MIRAS_front     = (0.30*0.1 + eps_mli*0.9);
alpha_MIRAS_back    = alpha_mli;
eps_MIRAS_back      = eps_mli;

% PROTEUS (??? same material as MIRAS gold coating)
A_bus               = 1;
alpha_bus           = alpha_mli; 
eps_bus             = eps_mli;

% Hub
A_hub               = 0.78;
A_hex               = 1.098;
alpha_hub           = alpha_mli;
eps_hub             = eps_mli;

% TOTAL QUANTITIES
A_tot       = A_arms*2 + A_bus*5 + A_hub*6 + A_hex;
r_sphere    = sqrt(A_tot/(4*pi));
A_cross_sun = pi*r_sphere^2;

alpha_SC = ((A_arms+A_hex)*alpha_MIRAS_front + A_arms*alpha_MIRAS_back + ...
            A_bus*5*alpha_bus + A_hub*6*alpha_hub) / A_tot;

eps_SC = ((A_arms+A_hex)*eps_MIRAS_front + A_arms*eps_MIRAS_back + ...
            A_bus*5*eps_bus + A_hub*6*eps_hub) / A_tot;

fprintf('\n========== SPACECRAFT PROPERTIES ==================\n')
fprintf(' - Averaged Absorbtivity   = %.3f \n', alpha_SC)
fprintf(' - Averaged Emissivity     = %.3f \n', eps_SC)
fprintf(' - Total Surface Area      = %.2f m^2\n', A_tot)

%% HEAT FLUX CALCULATIONS

% SUN
q_0         = 1367.5;
q_sun       = q_0;

% ALBEDO
a           = 0.4;
R_E         = 6371;                     % Earth radius [km]
h           = 765;                      % Orbit altitude [km]
r           = R_E + h;                  % Orbital radius [km]
q_alb       = q_sun * a * (R_E/r)^2;

% INFRARED
sigma       = 5.67e-8;
eps_E       = 0.8;
T_E         = 255;
q_ir        = sigma * eps_E * T_E^4 * (R_E/r)^2;

fprintf('\n========== HEAT FLUXES ============================\n')
fprintf(' - Solar                   = %.2f W/m^2\n', q_sun)
fprintf(' - Albedo                  = %.2f W/m^2\n', q_alb)
fprintf(' - Infrafed                = %.2f W/m^2\n', q_ir)

%% HEAT POWER CALCULATIONS

% SUN
Q_sun       = A_cross_sun * alpha_SC * q_sun;

% ALBEDO
K           = 1;                              %[TODO]
Q_alb       = A_tot * alpha_SC * K * q_alb;

% INFRARED
Q_ir        = A_tot * eps_SC * q_ir;

% INTERNAL
Q_int_mode  = [653.0-100; 761.0-190; 672.0-100; 710.2-120; 250.6-120];
Q_int_H     = max(Q_int_mode);
Q_int_C     = min(Q_int_mode);

fprintf('\n========== HEAT TRANSFER RATES ====================\n')
fprintf(' - Solar                   = %.2f W\n', Q_sun)
fprintf(' - Albedo                  = %.2f W\n', Q_alb)
fprintf(' - Infrafed                = %.2f W\n', Q_ir)
fprintf(' - Internal H              = %.2f W\n', Q_int_H)
fprintf(' - Internal C              = %.2f W\n', Q_int_C)

%% NON-REGULATED TEMPERATURES

% EMITTED
T_space     = 3;
% Q_emi     = sigma * eps_SC * A_tot * (T_SC^4 - T_space^4);

T_SC_hot    = ((Q_sun + Q_ir + Q_alb + Q_int_H)/(sigma*eps_SC*A_tot) + T_space^4)^0.25;
T_SC_cold   = ((Q_ir + Q_int_C)/(sigma*eps_SC*A_tot) + T_space^4)^0.25;

fprintf('\n========== NON-REG TEMPERATURES ===================\n')
fprintf(' - T_hot                   = %.1f °C\n', T_SC_hot-273.15)
fprintf(' - T_cold                  = %.1f °C\n', T_SC_cold-273.15)

%% RADIATORS
frac        = 0.5;
A_rad_max   = (A_arms + A_bus*4 + A_hub*2)*frac;
RAD_frac    = A_rad_max / A_tot;
eps_rad     = 0.85;
alpha_rad   = 0.10;

alpha_SC_rad = (...
    ( ...
    (A_arms+A_hex)*alpha_MIRAS_front + A_arms*(1-frac)*alpha_MIRAS_back + ...
    A_bus*2*alpha_bus + A_bus*3*(1-frac)*alpha_bus + ...
    A_hub*3*alpha_hub + A_hub*3*(1-frac)*alpha_hub ) + ...
    A_rad_max * alpha_rad) ...
    / A_tot;
eps_SC_rad = ( ...
    ( ...
    (A_arms+A_hex)*eps_MIRAS_front + A_arms*(1-frac)*eps_MIRAS_back + ...
    A_bus*2*eps_bus + A_bus*3*(1-frac)*eps_bus + ...
    A_hub*3*eps_hub + A_hub*3*(1-frac)*eps_hub ) + ...
    A_rad_max * eps_rad) ...
    / A_tot;

fprintf('\n========== with RADIATOR SC PROPERTIES ============\n')
fprintf(' - Max Radiator Area       = %.2f m^2 \n', A_rad_max)
fprintf(' - Radiator Fraction       = %.2f %% \n', RAD_frac*100)
fprintf(' - Radiator Absorptivity   = %.3f \n', alpha_rad)
fprintf(' - Radiator Emissivity     = %.3f \n', eps_rad)
fprintf(' - Averaged Absorptivity   = %.3f \n', alpha_SC_rad)
fprintf(' - Averaged Emissivity     = %.3f \n', eps_SC_rad)

%% RADIATOR HOT TEMPERATURE SOLUTION

Q_sun       = A_cross_sun * alpha_SC_rad * q_sun;
Q_alb       = A_tot * alpha_SC_rad * K * q_alb;
Q_ir        = A_tot * eps_SC_rad * q_ir;

T_SC_hot_rad = ((Q_sun + Q_ir + Q_alb + Q_int_H)...
    / (sigma*eps_SC_rad*A_tot) + T_space^4)^0.25;

fprintf('\n========== RADIATOR TEMP SOLUTION =================\n')
fprintf(' - T_hot_rad               = %.1f °C\n', T_SC_hot_rad-273.15)

%% HEATERS COLD TEMPERATURE SOLUTION

Q_emit = sigma*eps_SC_rad*A_tot * ((T_min_SC + 15)^4 - T_space^4);
Q_heat_req = (Q_emit - Q_ir - Q_int_C) * 1.25;

fprintf('\n======== HEATER T_min = %.1f °C SOLUTION =========\n', T_min_SC-273.15)
fprintf(' - Q Heaters + 25%% margin  = %.1f W\n', Q_heat_req)

%% TIME-VARIANT CASE

T_orb           = 100*60;
tspan           = linspace(0, 30*T_orb, 10000);
y0              = 273.15+10;

[t, y] = ode45( ...
    @(t, T) scThermalODE(t, T, Q_sun, Q_alb, Q_ir, Q_int_mode, eps_SC_rad), ...
    tspan, y0);

plot(t/T_orb, y-273.15, 'LineWidth', 2.5)
grid on

ax = gca;
ax.FontSize = 16; 

title('Time-variant thermal simulation', 'FontSize', 20)
xlabel('Orbits [-]', 'FontSize', 16); 
ylabel('Temperature [°C]', 'FontSize', 16)

yline(T_min_SC-273.15, 'r', 'LineWidth', 1.5);
yline(T_max_SC-273.15, 'r', 'LineWidth', 1.5);
ylim([T_min_SC-50 T_max_SC+10]-273.15)

legend('Temperature', 'Min. Temp.', 'Max. Temp.', 'Location', 'southwest', 'FontSize', 14)

fprintf('\n=========== TIME-VARIANT HEATER SOLUTION ==========\n')
fprintf(' - Min Temp. Reached       = %.1f °C\n', min(y)-273.15)
fprintf(' - Max Temp. Reached       = %.1f °C\n', max(y)-273.15)
