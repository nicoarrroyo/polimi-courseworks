function [V_in_list, V_out_list, dv_array, tof_array] = ...
    deep_space_injection(R_dep_list, V_dep_list, R_arr_list, ...
    dep_times, arr_times, steps, full_lambert, dv_lim)

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

V_in_list = NaN(steps, steps, 3);
V_out_list = NaN(steps, steps, 3);
dv_array = NaN(steps, steps, 3);
tof_array = NaN(steps, steps, 3);

for i = 1:steps
    temp_t1 = dep_times(i) * 24 * 3600;
    temp_R1 = R_dep_list(i, :);
    
    for j = 1:steps
        % check for arrival being after departure
        temp_t2 = arr_times(j) * 24 * 3600;
        temp_tof = temp_t2 - temp_t1;
        if temp_tof < 0
            continue
        end

        % check for valid lambert arc
        temp_R2 = R_arr_list(j, :);
        [~, ~, ~, temp_ERROR, temp_V1, temp_V2, ~, ~] = ...
            lambertMR(temp_R1, temp_R2, temp_tof, mu_sun, full_lambert, 0, 0, 0);
        if temp_ERROR ~= 0
            continue
        end

        % check for reasonable delta-v
        temp_dv = temp_V1 - V_dep_list(i, :);
        if norm(temp_dv) > dv_lim
            continue
        end

        V_in_list(i, j, :) = temp_V1;
        V_out_list(i, j, :) = temp_V2;
        dv_array(i, j, :) = temp_dv;
        tof_array(i, j) = temp_tof;
    end
end
