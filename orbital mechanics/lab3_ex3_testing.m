%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
labs_d = cd(1:backs(end)); addpath([labs_d '\student_functions']); 
addpath([labs_d '\lib']); addpath([labs_d '\lib' '\timeConversion']);
clear; close all; clc;

% Exercise 3: Mars Express
% Design an interplanetary transfer with minimum 𝚫𝒗𝐭𝐨𝐭 between Earth and 
% Mars, under the following mission requirements:
% time_depart_early = 2003 April 1; time_depart_late = 2003 August 1;
% time_arrive_early = 2003 September 1; time_arrive_late = 2004 March 1;
% No need to consider planetary insertion date.

%% 1. Implement a function to compute the 𝚫𝒗𝐭𝐨𝐭(t1, t2)
% see transfer_cost function in student_functions directory

%% 2. Evaluate Δ𝑣tot for a grid of departure and arrival times covering ...
% the time windows provided.

% Departure: 2003 April 1 - 2003 August 1
date_dep_min = [2003, 4, 1, 0, 0, 0];
date_dep_max = [2003, 8, 1, 0, 0, 0];

% Arrival: 2003 September 1 - 2004 March 1
date_arr_min = [2003, 9, 1, 0, 0, 0];
date_arr_max = [2004, 3, 1, 0, 0, 0];

date_dep_mjd2000 = date2mjd2000([2003, 06, 07, 22, 27, 34.14]);
date_arr_mjd2000 = date2mjd2000([2003, 12, 28, 14, 26, 08.25]);

earth_id = 3; mars_id = 4;

dvtot = transfer_cost(date_dep_mjd2000, date_arr_mjd2000, earth_id, mars_id);
tof_days = date_arr_mjd2000 - date_dep_mjd2000;



mu_sun = astroConstants(4);
[rE, vE] = get_planet_state(date_dep_mjd2000, earth_id, mu_sun);
[rM, vM] = get_planet_state(date_arr_mjd2000, mars_id, mu_sun);


% early_depart_mjd2000 = date2mjd2000([2003, 04, 01, 12, 00, 00]);
% early_arrive_mjd2000 = date2mjd2000([2003, 09, 01, 12, 00, 00]);
% 
% late_depart_mjd2000 = date2mjd2000([2003, 08, 01, 12, 00, 00]);
% late_arrive_mjd2000 = date2mjd2000([2004, 03, 01, 12, 00, 00]);
% 
% steps = 100;
% depart_times = linspace(early_depart_mjd2000, late_depart_mjd2000, steps);
% arrive_times = linspace(early_arrive_mjd2000, late_arrive_mjd2000, steps);
% 
% dvtot = zeros(length(depart_times), length(arrive_times));
% tof_days = zeros(length(depart_times), length(arrive_times));
% 
% for i = 1:length(depart_times)
%     for j = 1:length(arrive_times)
%         dvtot(i, j) = transfer_cost(...
%             depart_times(i), arrive_times(j), 3, 4);
% 
%         tof_days(i, j) = arrive_times(j) - depart_times(i);
%     end
% end
% 
% %% 3. Draw the porkchop plot of the Mars Express Mission
% figure("Name", "Mars Express Porkchop Plot"); hold on;
% 
% %v_levels = linspace(3, 10, 50);
% contourf(depart_times, arrive_times, dvtot');%, v_levels);
% 
% %tof_levels = 100:50:400;
% %contour(depart_times, arrive_times, tof_days, tof_levels, "ShowText", "on");
% 
% colorbar;
% xlabel("Departure Time (MJD2000)");
% ylabel("Arrival Time (MJD2000)");
% title("Porkchop Plot: Δv_{tot} for Mars Express Mission");
% hold off;
