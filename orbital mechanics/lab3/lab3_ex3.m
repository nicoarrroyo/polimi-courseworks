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
early_depart_mjd2000 = date2mjd2000([2003, 04, 01, 12, 00, 00]);
early_arrive_mjd2000 = date2mjd2000([2003, 09, 01, 12, 00, 00]);

late_depart_mjd2000 = date2mjd2000([2003, 08, 01, 12, 00, 00]);
late_arrive_mjd2000 = date2mjd2000([2004, 03, 01, 12, 00, 00]);

steps = 10;
depart_times = linspace(early_depart_mjd2000, late_depart_mjd2000, steps);
arrive_times = linspace(early_arrive_mjd2000, late_arrive_mjd2000, steps);

dvtot = zeros(length(depart_times), length(arrive_times));

tic
for i = 1:length(depart_times)
    for j = 1:length(arrive_times)
        dvtot(i, j) = transfer_cost(...
            depart_times(i), arrive_times(j), 3, 4);
    end
    fprintf("completed %.0f / %.0f departure calculations\n", ...
        i, length(depart_times))
end
toc

%% 3. Draw the porkchop plot of the Mars Express Mission
