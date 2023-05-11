function [T] = HomogeneousTransformOZ(oRelative, zAxisRelative)
% HOMOGENEOUSTRANSFORM - Finds the relative joint frame (from its parent)
% with the x axis along the common normal from parent z (0,0,1) to this z

% Inputs:
%   oRelative       - Frame origin, in parent coordinates
%   zAxisRelative   - Direction vector of z axis, in parent coordinates.
%                     Must be non-zero, but allowed to be non-unit length.
%
% Outputs:
%   T  - homogeneous transformation matrix.

% Authors: 
% Daniel Feshbach <feshbach@seas.upenn.edu>
% Last Edited 5/11/2023
%
% Copyright (C) 2023 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


zHat = zAxisRelative / norm(zAxisRelative);

% find common normal from parent
if zHat(1)~=0 || zHat(2)~=0
    commonNormal = cross([0,0,1], zHat);
else %zHat is parallel to parent Z, so check where origin is
    if oRelative(1)==0 && oRelative(2)==0 % Z axes coincident
        commonNormal = [1,0,0]; % no common normal, arbitrarily re-use parent X
    else
        commonNormal = [oRelative(1), oRelative(2), 0];
    end
end
xHat = commonNormal / norm(commonNormal);
yHat = cross(zHat, xHat);

T = zeros(4);
T(1:3,1) = xHat;
T(1:3,2) = yHat;
T(1:3,3) = zHat;
T(1:3,4) = oRelative;
T(4,4) = 1;

end