clc, clearvars, close all

%% TEMPERATURE RANGE

% Needs the real operating values of all the onboard components
% Tmax is the minimum of the maximum temperatures
% Tmin is the maximum of the minumum temperatures

T_max_SC            = 273.15 + 40;
T_min_SC            = 273.15 - 20;

%% AREA AND ALPHA/EPS CALCULATIONS (TODO)

% A_tot calculation needed ---> turn into equivalent sphere
% Different absorbptivity and emissivity for the SC surface
% Weighted average needed

% Solar Arrays (front: silicon solar cells | back: conductive/painted thermal finish)
A_sa                = 9.6;
alpha_sa_front      = 0.80;
eps_sa_front        = 0.82;
alpha_sa_back       = 0.90;
eps_sa_back         = 0.70;


A_tot       = A_sa*2;
A_cross_sun = A_sa;

alpha_SC = (A_sa*alpha_sa_front + A_sa*alpha_sa_back ) / A_tot;
eps_SC = (A_sa*eps_sa_front + A_sa*eps_sa_back ) / A_tot;

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
Q_sun       = A_cross_sun * alpha_sa_front * q_sun;

% ALBEDO
K           = 1;
Q_alb       = A_tot * alpha_SC * K * q_alb;

% INFRARED
Q_ir        = A_tot * eps_SC * q_ir;

% EMITTED
T_space     = 3;
% Q_emi      = sigma * eps_SC * A_tot * (T_SC^4 - T_space^4);

fprintf('\n========== HEAT TRANSFER RATES ====================\n')
fprintf(' - Solar                   = %.2f W\n', Q_sun)
fprintf(' - Albedo                  = %.2f W\n', Q_alb)
fprintf(' - Infrafed                = %.2f W\n', Q_ir)

%% NON-REGULATED TEMPERATURES

T_SC_hot    = ((Q_sun + Q_ir + Q_alb)/(sigma*eps_SC*A_tot) + T_space^4)^0.25;
T_SC_cold   = ((Q_ir)/(sigma*eps_SC*A_tot) + T_space^4)^0.25;

fprintf('\n========== NON-REG TEMPERATURES ===================\n')
fprintf(' - T_hot                   = %.1f °C\n', T_SC_hot-273.15)
fprintf(' - T_cold                  = %.1f °C\n', T_SC_cold-273.15)

%% RADIATORS
frac        = 0.60;
A_rad_max   = (A_sa)*frac;
RAD_frac    = A_rad_max / A_tot;
eps_rad     = 0.85;
alpha_rad   = 0.1;

alpha_SC_rad = (A_sa*alpha_sa_front + A_sa*alpha_sa_back*(1-frac)  + ...
            A_rad_max * alpha_rad) / A_tot;
eps_SC_rad = (A_sa*eps_sa_front + A_sa*eps_sa_back*(1-frac) + ...
            A_rad_max * eps_rad) / A_tot;

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

T_SC_hot_rad = ((Q_sun + Q_ir + Q_alb)...
    / (sigma*eps_SC_rad*A_tot) + T_space^4)^0.25;

fprintf('\n========== RADIATOR TEMP SOLUTION =================\n')
fprintf(' - T_hot_rad               = %.1f °C\n', T_SC_hot_rad-273.15)

%% HEATERS COLD TEMPERATURE SOLUTION

Q_emit = sigma*eps_SC_rad*A_tot * ((T_min_SC + 15)^4 - T_space^4);
Q_heat_req = (Q_emit - Q_ir) * 1.25;

fprintf('\n======== HEATER T_min = %.1f °C SOLUTION =========\n', T_min_SC-273.15)
fprintf(' - Q Heaters + 25%% margin  = %.1f W\n', Q_heat_req)

%% TIME-VARIANT CASE

T_orb           = 100*60;
tspan           = linspace(0, 30*T_orb, 10000);
y0              = 273.15+10;

[t, y] = ode45( ...
    @(t, T) scThermalODE(t, T, Q_sun, Q_alb, Q_ir, Q_int_mode, eps_SC_rad), ...
    tspan, y0);

plot(t/T_orb, y-273.15)
grid on
title('Time-variant thermal simulation')
xlabel('Orbits [-]'); ylabel('Temperature [°C]')
yline(T_min_SC-273.15, color='r'); yline(T_max_SC-273.15, color='r');
ylim([T_min_SC-50 T_max_SC+10]-273.15)
legend('Temperature', 'Min. Temp.', 'Max. Temp.')

fprintf('\n=========== TIME-VARIANT HEATER SOLUTION ==========\n')
fprintf(' - Min. Temp. Reached      = %.1f °C\n', min(y)-273.15)
