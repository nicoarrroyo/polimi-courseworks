function v_rotated = rodrigues(v, u, delta)
%   rodrigues.m Rodrigues' Formula to rotate a vector v by an angle 
%   delta around unit vector u (counter-clockwise)
%
%   SYNTAX
%   ------
%   v_rotated = functionName(u, v, delta)
%
%   INPUT ARGUMENTS
%   ---------------
%   v               - (3x1 matrix) vector to rotate.
%
%   u               - (3x1 matrix) vector around which to rotate.
%
%   delta           - (scalar) Rotation angle [rad].
%
%   OUTPUT ARGUMENTS
%   ----------------
%   v_rotated       - (3x1 matrix) rotated velocity matrix.
%
%   NOTES
%   -----
% 	Usage: angle delta must be in radians
% 	Dependencies: none
%
%   Author: Nicolas Renato Arroyo
%   Date: 2025 - 12 - 12
%   Updated: 2025 - 12 - 12

v_rotated = ...
    v * cos(delta) ...
    + cross(u, v) * sin(delta) ...
    + u * dot(u, v) * (1 - cos(delta));


end
