function dTdt = scThermalODE(t, T, Q_sun, Q_alb, Q_ir, Q_int_mode, eps_SC_rad)

C_SC = 5.6*10^5;

% ECLIPSE CHECK
T_orb = 100*60;
ecl_per = 17.5*60;
sun_per = T_orb - ecl_per;

if mod(t,T_orb) <= sun_per
    ecl_flag = 0;
else
    ecl_flag = 1;
end

% HEAT BALANCE
sigma = 5.67e-8;
A_tot = 18.94;
T_space = 3;
Q_heaters_max = 0;

if ecl_flag
    Q_in = Q_ir + min(Q_int_mode) + Q_heaters_max;
else
    Q_in = Q_sun + Q_alb + Q_ir + max(Q_int_mode);
end

Q_out = sigma*eps_SC_rad*A_tot*(T^4 - T_space^4);
dTdt = (Q_in - Q_out)/C_SC;
