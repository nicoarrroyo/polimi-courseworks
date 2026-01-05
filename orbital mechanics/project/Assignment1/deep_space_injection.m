function [V_out_list, V_in_list, dv_array, tof_array] = ...
    deep_space_injection(R_dep_list, V_dep_list, R_arr_list, V_arr_list, ...
    dep_times, arr_times, steps, dv_lim)

% FUNCTIONNAME Brief one-line description of what the function does.
%
%   SYNTAX
%   ------
%   output = functionName(input1, input2)
%   output = functionName(input1, input2, "OptionName", OptionValue)
%
%   DESCRIPTION
%   -----------
%   More detailed description of the function's purpose and behavior.
%   Explain what the function does, any important notes about its operation,
%   and any special cases or limitations.
%
%   INPUT ARGUMENTS
%   ---------------
%   input1			- (dataType) Description of input1.
%            			Default: [default value if applicable]
%
%   input2   		- (dataType) Description of input2.
%                  		Default: [default value if applicable]
%
%   OUTPUT ARGUMENTS
%   ----------------
%   output      	- (dataType) Description of output.
%
%   OPTIONAL PARAMETERS
%   -------------------
%   "OptionName" 	- (dataType) Description of the option.
%                  		Default: [default value]
%
%   EXAMPLES
%   --------
%   Example 1: Basic usage
%       output = functionName(input1, input2);
%
%   Example 2: With optional parameters
%       output = functionName(input1, input2, "OptionName", value);
%
%   NOTES
%   -----
%	Performance: 
% 	Numerical Stability: 
% 	Usage: 
% 	Dependencies: 
%
%   Author:  
%   Date:    
%   Updated: 
%
%   See also RELATEDFUNCTION1, RELATEDFUNCTION2

mu_sun = astroConstants(4);

V_out_list = NaN(steps, steps, 3);
V_in_list = NaN(steps, steps, 3);
dv_array = NaN(steps, steps, 3);
tof_array = NaN(steps, steps, 3);

for i = 1:steps
    t1 = dep_times(i) * 24 * 3600;
    R1 = R_dep_list(i, :);
    
    for j = 1:steps
        t2 = arr_times(j) * 24 * 3600;
        tof = t2 - t1;

        % check for arrival being after departure
        if tof <= (30*24*3600) % minimum time for tof
            continue
        end

        R2 = R_arr_list(j, :);
        [~, ~, ~, ERROR, V1, V2, ~, ~] = ...
            lambertMR(R1, R2, tof, mu_sun, 0, 0, 0, 0);
        
        % check for valid lambert arc
        if ERROR ~= 0
            continue
        end

        dv = V1 - V_dep_list(i, :);

        % check for reasonable delta-v
        if norm(dv) > dv_lim
            continue
        end

        V_out_list(i, j, :) = V1;
        V_in_list(i, j, :) = V2;
        dv_array(i, j, :) = dv;
        tof_array(i, j) = tof;
    end
end
