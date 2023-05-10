function [DataStruct] = SphereMinimization(DataStruct, index, r, colorinput, plotoption)
% SPHEREMINIMIZATION - Optimization function which reassigns the
% centerpoint of a given bounding sphere, provided that it is allowed to
% move along its z axis and must be a specified distance from previous
% bounding spheres.

% Inputs:
%   DataStruct - array which contains information about the origins of each
%                of the joints as well as the radius associated with each
%                joint. We will be starting with the end joint and working
%                to joint 0.
%   index      - value which allows us to keep track of what joint we are 
%                adjusting. index will be counted in an external loop.
%   r          - the radius unadjusted by any joint parameters (pure r).
%   colorinput - string input which determines the color of the plotted
%                sphere, if so desired.
%   plotoption - string input which dictates plotting.

% Outputs:
%   DataStruct - updated array which contains information about the origins
%                of each of the joints as well as the radius associated
%                with each joint. 

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Daniel Feshbach <feshbach@seas.upenn.edu>
% Last Edited 4/18/2023
%
% Copyright (C) 2023 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% Calculate required distances to make conditionals easier to write
dist = 4*r + DataStruct(index).rnew + DataStruct(index-1).rs;
% x_dist = abs(DataStruct(index).oinew(1) - DataStruct(index-1).oi(1));
% y_dist = abs(DataStruct(index).oinew(2) - DataStruct(index-1).oi(2));
% z_dist = abs(DataStruct(index).oi(3) - DataStruct(index-1).oi(3));

zzz = -DataStruct(index-1).zglobal;
% zzz = [4,1,1];
% Choose delta such that norm(DataStruct(index-1).oinew -
% DataStruct(index).oinew) is minimized but also so that
% norm(DataStruct(index-1).oinew - DataStruct(index).oilarge) >= dist

% DataStruct(index-1).oinew = DataStruct(index-1).oi + zzz*delta;

% Statement for minimization
obj = @(delta) norm(DataStruct(index-1).oi + zzz*delta - DataStruct(index).oinew);

% Initial Guess
delta0 = 1;

% Inequality requirements:
% norm(DataStruct(index-1).oi + zzz*delta - DataStruct(index).oi) >= dist
% Nonlinear

    % Nonlinearity constraint
    function [c, ceq] = nlcon(delta)
        c = dist - norm(DataStruct(index-1).oi + zzz*delta - DataStruct(index).oilarge);
        ceq = [];
    end


if isnan(DataStruct(index-1).delta_z)   
    % Express other parameters as empty cells (not used)
    A = [];
    B = [];
    Aeq = [];
    Beq = [];
    LB = [];
    UB = [];
    nonlincon = @nlcon;
    
    % Output optimal value for delta
    delta_opt = fmincon(obj, delta0, A, B, Aeq, Beq, LB, UB, nonlincon);
    
    % negative sign because zzz is defined as -DataStruct(index-1).zglobal;
    DataStruct(index-1).delta_z = -delta_opt;
end

% Given delta_opt, determine new value for DataStruct(index-1).oinew
DataStruct(index-1).oinew = DataStruct(index-1).oi - zzz * DataStruct(index-1).delta_z;

% hold on
% quiver3(DataStruct(index-1).oinew(1), DataStruct(index-1).oinew(2), ...
%     DataStruct(index-1).oinew(3), zzz(1), zzz(2), zzz(3), ...
%     'AutoScaleFactor', 3, 'Linewidth', 3);

% dist - norm(x - DataStruct(index).oi)
% 
% p_oi = DataStruct(index-1).oi + zzz*delta;
% p_oinew = DataStruct(index-1).oinew;
% prob = optimproblem;
% prob.Objective = norm(p_oinew - p_oi);
% prob.Constraints.cons1 = norm(p_oinew - p_oi) > dist;
% 
% sol = solve(prob);
% 
% DataStruct(index-1).oinew = 

% Solve delta such that norm(p_oinew - p_oi) is min, but  norm(p_oinew - p_oi) > dist
% Optimization

% if x_dist >= dist || y_dist >= dist
%     
%     DataStruct(index-1).oinew = DataStruct(index-1).oi;
%     
% elseif sqrt(x_dist^2 + y_dist^2) >= dist
%         
%     DataStruct(index-1).oinew = DataStruct(index-1).oi;
%         
% else
%     
%     if x_dist > y_dist
%         
%         % Use x term if it is greater, as it will result in a smaller z
%         % term required for offset (beneficial)
%         DataStruct(index-1).oinew(3) = DataStruct(index).oinew(3) - sqrt(dist^2 - x_dist^2);
%         DataStruct(index-1).oinew(1) = DataStruct(index-1).oi(1);
%         DataStruct(index-1).oinew(2) = DataStruct(index-1).oi(2);
%           
%     else
%         
%         % Otherwise use y term for same reason
%         DataStruct(index-1).oinew(3) = DataStruct(index).oinew(3) - sqrt(dist^2 - y_dist^2);
%         DataStruct(index-1).oinew(1) = DataStruct(index-1).oi(1);
%         DataStruct(index-1).oinew(2) = DataStruct(index-1).oi(2);
%         
%     end
%         
% end

for i = size(DataStruct, 2):-1:index-1
    
    [DataStruct(i).demo] = SphericalSampling(DataStruct(i).oinew, ...
        DataStruct(i).rs, colorinput, plotoption);
    
    % If on the last index, the cumulative sphere points is just the final
    % sphere itself
    if i == size(DataStruct, 2)
        
        DataStruct(i).democumul = DataStruct(i).demo;  
    else
        
        DataStruct(i).democumul = [DataStruct(i).demo; DataStruct(i+1).democumul];
    end
    
end

% Find new sphere center and store to index-1.optimized
[R, C, Xb] = ExactMinBoundSphere3D(DataStruct(index-1).democumul);

DataStruct(index-1).Xb = Xb;
DataStruct(index-1).oilarge = C;
DataStruct(index-1).rnew = R;

% Output new plot of cumulative sphere
[DataStruct(index-1).spheresnet] = SphericalSampling(DataStruct(index-1).oilarge, ...
    DataStruct(index-1).rnew, colorinput, plotoption);

      
end