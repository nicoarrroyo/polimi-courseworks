function dy = ode_2bp_j2( ~, y, mu, J2, R_e )

r = y(1:3);
v = y(4:6);

x = r(1);
y = r(2);
z = r(3);

r_norm = norm(r);

a_j2_x = (x / r_norm) * ((5 * z^2 / r_norm^2) - 1);
a_j2_y = (y / r_norm) * ((5 * z^2 / r_norm^2) - 1);
a_j2_z = (z / r_norm) * ((5 * z^2 / r_norm^2) - 3);

a_j2 = (3/2) * (J2 * mu * R_e^2 / r_norm^4) * [a_j2_x; a_j2_y; a_j2_z;];

dy = [  v 
        (-mu/r_norm^3)*r + a_j2 ];

end