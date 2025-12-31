%% configure paths
cd = fileparts(mfilename("fullpath")); backs = strfind(cd, "\"); 
proj_d = cd(1:backs(end)); addpath([proj_d '\student_functions']); 
addpath([proj_d '\lib']); addpath([proj_d '\lib' '\timeConversion']);
clear; close all; clc;

%% === Assignment 1 ===

%% --- Mission Requirements ---
% --- Travel Window ---
travel_window.start_date = [2030, 1, 1, 0, 0, 0];
travel_window.end_date = [2060, 1, 1, 0, 0, 0];
travel_window.start_mjd2k = date2mjd2000(travel_window.start_date);
travel_window.end_mjd2k = date2mjd2000(travel_window.end_date);

travel_window.tof = travel_window.start_date - travel_window.end_date;

% --- Departure Planet ---
planet_dep.name = "Mercury";
planet_dep.mu = astroConstants(11);

% --- Flyby Planet ---
planet_fb.name = "Earth";
planet_fb.mu = astroConstants(13);

% --- Arrival Asteroid ---
asteroid_arr.id = 316801;
asteroid_arr.name = "N." + asteroid_arr.id;
[...
    asteroid_arr.kep1, ...
    asteroid_arr.mass, ...
    asteroid_arr.M1...
    ] = ephAsteroids(travel_window.start_mjd2k, asteroid_arr.id);
[...
    asteroid_arr.kep2, ...
    ~, ...
    asteroid_arr.M2...
    ] = ephAsteroids(travel_window.end_mjd2k, asteroid_arr.id);
