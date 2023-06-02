function [theta] = SignedAngle(a, b, n)
% SIGNEDANGLE - Finding the signed angle between two vectors "a" and "b"
% with respect to vector "a", around the normal vector "n". Outputs an
% angle "theta".

% Inputs:
%   a       - vector "a". See function description for details.
%   b       - vector "b". See function description for details.
%   n       - normal vector for utilization in cross product.

% Outputs:
%   theta   - output signed angle, in radians.

% Authors: 
% Wei-Hsi Chen <weicc@seas.upenn.edu>
% Daniel Feshbach <feshbach@seas.upenn.edu>
% Last edited 5/22/2023
%
% Copyright (C) 2023 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% Make sure all vectors are normalised
a = a/norm(a);
b = b/norm(b);
if norm(n)~=0
    n = n/norm(n);
end


% find the signed angle
theta = atan2(dot(cross(a,b),n),dot(a,b));

end

