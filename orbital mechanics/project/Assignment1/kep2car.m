function [r0, v0] = kep2car( a, e, i, Omega, omega, TA, mu, angle )
    % ----- INPUTS -----
    % a: semi-major axis
    % e: eccentricity scalar
    % i: inclination
    % RAAN: right ascension of ascending node
    % w (small omega): argument of perigee
    % TA: true anomaly
    % ------------------
    
    % ----- OUTPUTS -----
    % r0: cartesian position vector
    % v0: cartesian velocity vector
    % -------------------

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
    
    %% initial constants
    if nargin < 7
        mu = 398600;
    end
    if nargin < 8
        angle = "degrees";
    end
    % Convert angles to degrees if necessary
    if angle == "radians"
        i = rad2deg(i);
        Omega = rad2deg(Omega);
        omega = rad2deg(omega);
        TA = rad2deg(TA);
    end


    p = a * (1 - e^2);
    r = p / (1 + e * cosd(TA));
    h = sqrt(mu * p);
    
    %% keplerian to perifocal (eph)
    % radius
    r_eph = [r * cosd(TA); ...
             r * sind(TA); ...
             0;]; % perifocal radius vector
    
    % velocity
    vr = (mu / h) * e * sind(TA); % radial velocity component
    vt = (mu / h) * (1 + e * cosd(TA)); % transverse velocity component
    v_eph = [vr * cosd(TA) - vt * sind(TA); ...
             vr * sind(TA) + vt * cosd(TA); ...
             0;]; % perifocal velocity vector
    
    %% perifocal to cartesian
    % use a 3-1-3 rotation matrix. in order:
        % 1. rotate around Z by the argument of perigee (omega)
        % 2. rotate around X by the inclination (i)
        % 3. rotate around Z by the RAAN (Omega)
    R3_omega = [cosd(omega) -sind(omega) 0;
                sind(omega) cosd(omega) 0;
                0 0 1];
    
    R1_i = [1 0 0;
            0 cosd(i) -sind(i);
            0 sind(i) cosd(i)];
    
    R3_Omega = [cosd(Omega) -sind(Omega) 0;
                sind(Omega) cosd(Omega) 0;
                0 0 1];
    
    r0 = R3_Omega * R1_i * R3_omega * r_eph;
    v0 = R3_Omega * R1_i * R3_omega * v_eph;
end
