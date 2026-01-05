function [r, v] = get_planet_state(mjd2000, body_id, mu)
% GET_PLANET_STATE Get position and velocity of a planet
%
% Inputs:
%   mjd2000  - Modified Julian Date 2000 [days]
%   body_id  - Celestial body identifier
%   mu       - gravitational parameter of the body in question
%
% Outputs:
%   r - Position vector [km] (column vector)
%   v - Velocity vector [km s^-1] (column vector)

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

    [kep, ~] = uplanet(mjd2000, body_id);

    [r, v] = kep2car(...
        kep(1), ...  % Semi-major axis
        kep(2), ...  % Eccentricity
        kep(3), ...  % Inclination
        kep(4), ...  % RAAN
        kep(5), ...  % Argument of periapsis
        kep(6), ...  % True anomaly
        mu, ...
        "radians");
    
    r = r'; v = v'; % transpose for Lambert function
end