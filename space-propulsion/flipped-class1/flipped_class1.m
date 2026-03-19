clear; close all; clc;
% NOTES
% There is one important error in the data which you can and should find. 
% The medium and high curves are swapped for pbar2444
% Data is not cleaned. There is a spike at the beginning due to ignition. 
% This peak must be dealt with. 

%% --- load data --- %
cd("/home/nico/Documents/POLIMI/YEAR1/SEM2/SPACE PROPULSION/COURSEWORK")
load("02 - Non reacting flows - lecture slides-20260303/tracesbar1.mat");

%% --- apply bayern-chemie method --- %
pbar = pbar2438;
actionA = [0 0 0];
actionG = [0 0 0];
is_valid = 0;

for i = 1:3
    fprintf("\n============== COLUMN %.0f ==============\n", i)
    [max_p, max_p_idx] = max(pbar(:, i));
    found_up = 0;
    found_down = 0;
    for j = 1:length(pbar)
        if pbar(j, i) > (0.1*max_p) && found_up ~= 1
            found_up = 1;
            actionA(i) = pbar(j, i);
            
            fprintf("found action time A at %.0f\n", j)
            fprintf("maximum pressure %.2f\n", max_p)
            fprintf("upwards inflection point pressure %.2f\n", actionA(i))
        end
        if pbar(j, i) < (0.1*max_p) && found_down ~= 1 && found_up == 1
            found_down = 1;
            actionG(i) = pbar(j, i);
            
            fprintf("found action time G at %.0f\n", j)
            fprintf("maximum pressure %.2f\n", max_p)
            fprintf("downwards inflection point pressure %.2f\n", actionG(i))
        end
    end
    fprintf("======================================\n")
end

%% plot
figure()
plot(pbar(1:200, :), linewidth=1)
xlabel("time [ms]")
ylabel("pressure [bar]")
xline(actionA, "g")
xline(actionG, "r")
legend("low", "medium", "high", "A", "G")


%% --- determine vieille parameters --- %


%% --- calculate characteristic velocity --- %


%% --- build predictive motor model --- %


