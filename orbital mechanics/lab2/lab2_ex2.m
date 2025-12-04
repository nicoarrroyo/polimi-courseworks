%% configure paths
script_path = fileparts(mfilename("fullpath"));
backs = strfind(script_path, "\"); labs_dir = script_path(1:backs(end));
addpath([labs_dir '\student_functions']); addpath([labs_dir '\lib']);
clear; close all; clc;

clear; close all; clc;

%% constants / initial conditions
w_E = 15.04; % earth rotation velocity [deg hr^-1]
mu_E = astroConstants(13); % earth gravitational parameter [km^3 s^-2]
earth_img = imread("EarthTexture.jpg");

%% %%% generic orbit case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r0 = [-4578.219 -801.084 -7929.708]; % initial position vector [km]
v0 = [0.800 -6.037 1.385]; % initial velocity vector [km s^-1]
orbits = 3.25; % number of orbits to propogate for [-]
title = "Generic Orbit";
lab2_ground_track( r0, v0, orbits, title )

%% %%% molniya orbit case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r0 = [3108.128 -1040.299 -6090.022]; % initial position vector [km]
v0 = [5.743 8.055 1.555]; % initial velocity vector [km s^-1]
orbits = 30; % number of orbits to propogate for [-]
title = "Molniya Orbit";
lab2_ground_track( r0, v0, orbits, title )

%% %%% three circular LEO orbits case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r0 = [ % initial position vector [km]
    [5493.312 4609.436 0.000]
    [5493.312 3991.889 2304.718]
    [5493.312 -641.510 4564.578]
    ];
v0 = [ % initial velocity vector [km s^-1]
    [-4.792 5.711 0.000]
    [-4.792 4.946 2.856]
    [-4.792 -0.795 5.656]
    ];
k = [20; 29; 15;];
m = [2; 2; 1;];
i = [0; 30; 98;];
orbits = 5; % number of orbits to propogate for [-]
for i = 1:height(r0)
    title = "Circular LEO " + num2str(i);
    r_0 = r0(i, :);
    v_0 = v0(i, :);
    lab2_ground_track( r_0, v_0, orbits, title )
end
