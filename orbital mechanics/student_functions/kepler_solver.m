function E = kepler_solver( t, e, a, mu, t0, E0 )

n = sqrt(mu / a^3); % mean motion

M = n * (t - t0) + E0 - e * sin(E0);

kepler_eq = @(E) E - e * sin(E) - M;

E_guess = M + (e * sin(M)) / (1 - sin(M + e) + sin(M));

options = optimoptions("fsolve", "Display", "none", "TolFun", 1e-14);

E = fsolve(kepler_eq, E_guess, options);

end