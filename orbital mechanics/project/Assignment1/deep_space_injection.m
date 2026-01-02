function [V1_list, V2_list, dv_array, tof_array] = ...
    deep_space_injection(R_dep_list, V_dep_list, R_arr_list, ...
    V_arr_list, dep_times, arr_times, steps, full_lambert)

V1_list = zeros(steps, 3);
V2_list = zeros(steps, 3);
dv_array = NaN(steps, steps);
tof_array = NaN(steps, steps);

for i = 1:steps
    R1 = R_dep_list(i, :);
    V1 = V_dep_list(i, :);
    t1 = dep_times(i) * 24 * 3600;

    dv_row = NaN(1, steps);
    tof_row = NaN(1, steps);
    for j = 1:steps
        t2 = arr_times(j) * 24 * 3600;
        tof = t2 - t1;
        
        if tof > 0
            R2 = R_arr_list(j, :);

            [~, ~, ~, ERROR, V1_temp, V2_temp, ~, ~] = ...
                lambertMR(R1, R2, tof, mu_sun, 0, 0, 0, 0);
            if ERROR == 0
                tof_row(j) = t2 - t1;

                V1_list(j, :) = V1_temp;
                V2_list(j, :) = V2_temp;

                if full_lambert ~= 0 % if aiming for rendez-vous
                    V2 = V_arr_list(j, :);
                    dv_array(j) = norm(V1_temp - V1) + norm(V2_temp - V2);
                else % if aiming for only gravity assist
                    dv_array(j) = norm(V1_temp - V1);
                end
            end
        end
    end
    dv_array(i, :) = dv_row;
    tof_array(i, :) = tof_row / (24 * 3600);
end
