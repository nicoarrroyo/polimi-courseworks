clear; close all; clc;
% NOTES
% There is one important error in the data which you can and should find. 
% The medium and high curves are swapped for pbar2444
% Data is not cleaned. There is a spike at the beginning due to ignition. 
% This peak must be dealt with. 

%% --- load data --- %
cd("/home/nico/Documents/POLIMI/YEAR1/SEM2/SPACE PROPULSION/COURSEWORK")
load("02 - Non reacting flows - lecture slides-20260303/tracesbar1.mat");

% --- apply bayern-chemie method --- %
%% get action times A and D
pbar = pbar2439;
actionA_p = [NaN NaN NaN];
actionG_p = [NaN NaN NaN];
actionA_idx = [NaN NaN NaN];
actionG_idx = [NaN NaN NaN];

for i = 1:3
    fprintf("\n============== COLUMN %.0f ==============\n", i)
    [pmax, pmax_idx] = max(pbar(:, i));
    found_up = 0;
    found_down = 0;
    for j = 1:pmax_idx-1
        if pbar(pmax_idx-j, i) < 0.05*pmax && found_up ~= 1
            found_up = 1;
            actionA_p(i) = pbar(pmax_idx-j, i);
            actionA_idx(i) = pmax_idx-j;

            fprintf("action time A at %.0f ms\n", (pmax_idx-j))
            fprintf("pressure at A %.2f bar\n", actionA_p(i))
        end
    end
    for j = pmax_idx:length(pbar(:, i))
        if pbar(j, i) < 0.05*pmax && found_down ~= 1
            found_down = 1;
            actionG_p(i) = pbar(j, i);
            actionG_idx(i) = j;

            fprintf("action time G at %.0f ms\n", j)
            fprintf("pressure at G %.2f bar\n", actionG_p(i))
        end
    end
    fprintf("======================================\n")
end

%% get p_ref
p_ref = [NaN NaN NaN];

action_time = actionG_idx - actionA_idx;
for i = 1:3
    sum = 0;
    for j = actionA_idx(i):actionG_idx(i)-1
        sum = sum + 0.5*(pbar(j+1, i) + pbar(j, i));
    end
    p_ref(i) = sum / 2 / action_time(i);
end

%% get action times B and E
actionB_idx = [NaN NaN NaN];
actionE_idx = [NaN NaN NaN];

for i = 1:3
    fprintf("\n============== COLUMN %.0f ==============\n", i)
    [pmax, pmax_idx] = max(pbar(:, i));
    found_up = 0;
    found_down = 0;
    for j = 1:pmax_idx-1
        if pbar(pmax_idx-j, i) < p_ref(i) && found_up ~= 1
            found_up = 1;
            actionB_idx(i) = pmax_idx-j;

            fprintf("action time B at %.0f ms\n", (pmax_idx-j))
        end
    end
    for j = pmax_idx:length(pbar(:, i))
        if pbar(j, i) < p_ref(i) && found_down ~= 1
            found_down = 1;
            actionE_idx(i) = j;

            fprintf("action time E at %.0f ms\n", j)
        end
    end
    fprintf("======================================\n")
end

%% get p_eff
p_eff = [NaN NaN NaN];

burn_time = actionE_idx - actionB_idx;
for i = 1:3
    sum = 0;
    for j = actionB_idx(i):actionE_idx(i)-1
        sum = sum + 0.5*(pbar(j+1, i) + pbar(j, i));
    end
    p_eff(i) = sum / burn_time(i);
end

%% plot
% figure()
% plot(pbar(1:end, :), linewidth=1)
% xlabel("time [ms]")
% ylabel("pressure [bar]")
% xline(actionA_p(1), "g")
% xline(actionG_p(1), "r")
% legend("low", "medium", "high", "A-low", "G-low")

%% --- determine vieille parameters --- %


%% --- calculate characteristic velocity --- %


%% --- build predictive motor model --- %


