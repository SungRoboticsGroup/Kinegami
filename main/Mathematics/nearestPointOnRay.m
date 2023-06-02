function [nearestPoint] = nearestPointOnRay(p, rayStart, rayDirection)
% Inputs:
%   p             - a 3D point
%   rayStart      - a 3D point where the ray starts
%   rayDirection  - a non-zero 3D vector defining the direction of the ray

% Outputs:
%   nearestPoint  - the 3D point closest to pt on the line through
%                   pointOnLine in direction lineDirection


rayDirection = rayDirection / norm(rayDirection);
projectedDistanceAlongRay = dot((p - rayStart), rayDirection);
if projectedDistanceAlongRay <= 0
    nearestPoint = rayStart;
else
    nearestPoint = rayStart + projectedDistanceAlongRay*rayDirection;
end

end