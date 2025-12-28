function ds = eq_motion(t, s, acc_pert_fun, parameters)
acc_pert_vec = acc_pert_fun(t, s, parameters);

ds = ...

end