function porkchop_plot(plot_name, dv_array, dep_times, arr_times)

figure("Name", plot_name); hold on;

% --- plot dv ---
dv_min = min(min(dv_array)); dv_max = max(max(dv_array));
dv_step = (dv_max - dv_min) / 30;
v_levels = dv_min : dv_step : dv_max;
contour(dep_times, arr_times, dv_array', v_levels, "ShowText", "off");
clim([dv_min, dv_max]);

% --- plot the constant time of flight lines ---
% leg1.tof_levels = 50 : 50 : 200;
% contour(leg1.dep_times, leg1.arr_times, leg1.tof, leg1.tof_levels, "ShowText", "on");

% --- plot ---
colorbar; grid on;
xlabel("Departure Time (MJD2000)");
ylabel("Arrival Time (MJD2000)");
title("Δv_{tot}");
hold off;

end