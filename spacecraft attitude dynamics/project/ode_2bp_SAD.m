function dy = ode_2bp_SAD( ~, y, mu )

r = y(1:3);
v = y(4:6);

rnorm = norm(r);

dy = [  v 
        (-mu/rnorm^3)*r ];

end