function find_lowest_dv_mission(dv_array, dep_times, arr_times)

[min_val, min_idx_linear] = min(dv_array(:));
[row_idx, col_idx] = ind2sub(size(dv_array), min_idx_linear);

mjd_opt_dep_grid = dep_times(row_idx);
date_d = mjd20002date(mjd_opt_dep_grid);

mjd_opt_arr_grid = arr_times(col_idx);
date_a = mjd20002date(mjd_opt_arr_grid);

tof_opt_grid = mjd_opt_arr_grid - mjd_opt_dep_grid;

% --- results output ---
fprintf("=== GRID SEARCH RESULTS ===\n");
fprintf("Min Delta V: %.4f km s^-1\n", min_val);
fprintf("Departure:   MJD2000 %.2f\n", mjd_opt_dep_grid);
fprintf("             Date    %.0f %.0f %.0f %.0f %.0f %.0f\n", date_d);
fprintf("             ROW     %0.f\n", row_idx);
fprintf("Arrival:     MJD2000 %.2f\n", mjd_opt_arr_grid);
fprintf("             Date    %.0f %.0f %.0f %.0f %.0f %.0f\n", date_a);
fprintf("             COL     %0.f\n", col_idx);
fprintf("Total ToF (days):      %.2f\n", tof_opt_grid);
end