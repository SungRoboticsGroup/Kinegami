function nhat = unitNormalToBoth(u, v)
% UNITNORMALTOBOTH - Find a unit vector orthogonal to both u and v, 
%                    obeying the right-hand-rule if applicable (i.e., if
%                    they aren't colinear and neither is 0).

% Input:
% u, v      3D vectors

% Output
% nhat      A 3D unit column vector which dots to 0 with both u and v.
%           If u,v have nonzero cross product, nhat is the unit vector
%           in direction cross(u,v).

% Authors: 
% Daniel Feshbach <feshbach@seas.upenn.edu>
% Last edited 6/12/2023
%
% Copyright (C) 2023 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.

cp = cross(u,v);
if norm(cp) ~= 0
    nhat = cp / norm(cp);
    if isrow(nhat)
        nhat = nhat'; % want to return a column vector
    end
else
    % ensure both u and v are row vectors (needed to find nullspace)
    if ~isrow(u)
        u = u';
    end
    if ~isrow(v)
        v = v';
    end
    B = null([u;v]); % orthonormal basis for the null space
    b1 = B(:,1);     % the first basis vector 
    % (any nullspace basis vector will work, but we're not guaranteed 
    % more than one exists)
    nhat = b1;   
end

end




    % Old version from before I learned of the null function
    % % Input: 3d vectors u, v
    % % Output: a unit vector nhat normal to both u and v
    % %         works even if u and v are colinear
    % function nhat = unitNormalToBoth(u, v)
    %     % if either input is 0, we don't need 
    %     if norm(u) == 0
    %         u = [0,0,1];
    %     end
    %     if norm(v) == 0
    %         v = [0,1,0];
    %     end
    % 
    %     uxv = cross(u, v);
    %     if norm(uxv) == 0
    %         ux100 = cross(u, [1,0,0]);
    %         if norm(ux100) == 0
    %             % this means u,v are colinear with [1,0,0]
    %             % so they're normal to [0,1,0]
    %             nhat = [0,1,0];
    %         else
    %             nhat = ux100 / norm(ux100);
    %         end
    %     else
    %         nhat = uxv / norm(uxv);
    %     end
    % end