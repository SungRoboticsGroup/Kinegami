function [infostruct] = DubinsTube(r, n, Op, Od, infostruct, index, mirror, split)
% DUBINSTUBE - Generates a sequence of origami modules to connect the 
% distal frame of the previous joint and the proximal frame of the
% subsequent joint.

% Inputs:
%   r           - desired radius of folded origami linkage.
%   n           - number of sides of folded origami linkage.
%   Op          - the matrix which describes the proximal frame of the
%                 linkage in question. This is a 3x4 matrix so as to
%                 include the center point of the x, y, and z vectors.
%   Od          - the matrix which describes the distal frame of the
%                 linkage in question. This is a 3x4 matrix so as to
%                 include the center point of the x, y, and z vectors.
%   infostruct  - a data structure that includes all the information needed
%                 for the construction of the full schematic.
%   index       - used for indexing of outputs in infostruct. 
%   mirror      - setting required for creation of elbow schematic. See
%                 Origami_Elbow_CreasePattern.m for further detail.
%   split       - setting required for creation of elbow schematic. See
%                 Origami_Elbow_CreasePattern.m for further detail.

% Outputs:
%   infostruct  - updated data structure that includes all the information
%                 needed for the construction of the full schematic.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Wei-Hsi Chen <weicc@seas.upenn.edu>
% Last edited 7/11/2022
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% We know that Ou = {au, bu, cu, ou} - 3x4 matrix of vec/norm(vec)
% Thus,
ap = Op(:, 1);
ap = ap / norm(ap);

ad = Od(:, 1);
ad = ad / norm(ad);

bp = Op(:, 2);
bp = bp / norm(bp);

bd = Od(:, 2);
bd = bd / norm(bd);

% Run solution solver to compare four sets of solutions and find t with the
% shortest path. output this t and the corresponding theta1 and theta2
[Tdir, Tmag, theta1, theta2] = solveDubins3d(r, Od, Op);

tunit = Tdir/norm(Tdir);
t = Tmag * tunit;
infostruct(index).t = t;
infostruct(index).tunit = tunit;

% Determine values of wp and wd
wp = cross(ap, tunit);
wp = wp / norm(wp);

wd = cross(tunit, ad);
wd = wd / norm(wd);

% Define bm
[r_mat_p] = RotationalMatrix(wp, theta1);
bm = r_mat_p * bp;

% Determine values of phi1
% phi1 = atan2(norm(cross(bp, wp)), dot(bp, wp));
phi1 = SignedAngle(bp, wp, ap);

[r_mat_d] = RotationalMatrix(wd, theta2);

% Use results from rotational matrices, along with b_p and b_d, to
% determine the twist angle, alpha. 
bu = r_mat_d * bm;

% Use value of b_u to find alpha
% alpha = atan2(norm(cross(bu, bd)), dot(bu, bd));
alpha = SignedAngle(bu, bd, ad);

disp(alpha)

% Determine phi2
% phi2 = atan2(norm(cross(bm, wd)), dot(bm, wd)) - alpha;
phi2 = SignedAngle(bm, wd, tunit) - alpha;

% Elbow Fitting
[lengths, ls] = Origami_Elbow_Parameters(r, n, phi1, theta1, split);
[tuckangle] = TuckAngles(r, n, phi1, theta1, split);

infostruct(index).ls = ls;
infostruct(index).r = r;

[dataFoldA, m, lmax] = Origami_Elbow_CreasePattern(lengths, ls, n, infostruct(index).h1, ...
    infostruct(index).h2, theta1, mirror, split, tuckangle);

infostruct(index).m = m;
infostruct(index).lmax = lmax;
infostruct(index).n = n;
infostruct(index).type = dataFoldA;
infostruct(index).name = "Elbow";
infostruct(index).theta = theta1;
infostruct(index).dw = r*abs(tan(theta1 / 2));
infostruct(index).oc = Op(:, 4) + r*abs(tan(theta1 / 2))*ap;

% Twist Fitting
Tnorm = norm(t);

[x, l, ls] = Origami_Twist_Parameters(r, n, 0.2*Tnorm, alpha);

infostruct(index+1).ls = ls;
infostruct(index+1).r = r;

[dataFoldB, m, lmax] = Origami_Twist_CreasePattern(x, l, ls, n, infostruct(index+1).h1, ...
    infostruct(index+1).h2, alpha);

infostruct(index+1).m = m;
infostruct(index+1).lmax = lmax;
infostruct(index+1).n = n;
infostruct(index+1).type = dataFoldB;
infostruct(index+1).name = "Twist";
infostruct(index+1).alpha = alpha;
infostruct(index+1).h = 0.2*Tnorm;

% Tube
[ls] = Origami_Tube_Parameters(r, n);

infostruct(index+2).ls = ls;
infostruct(index+2).r = r;
 
% Outputs graphing for default tube
[dataFoldTube, m, lmax] = Origami_Tube_CreasePattern(n, ls, 0.8*Tnorm);

infostruct(index+2).m = m;
infostruct(index+2).lmax = lmax;
infostruct(index+2).n = n;
infostruct(index+2).type = dataFoldTube;
infostruct(index+2).name = "Tube";
infostruct(index+2).h = 0.8*Tnorm;

% Elbow Fitting pt. 2
[lengths, ls] = Origami_Elbow_Parameters(r, n, phi2, theta2, split);
[tuckangle] = TuckAngles(r, n, phi2, theta2, split);

infostruct(index+3).ls = ls;
infostruct(index+3).r = r;

[dataFoldD, m, lmax] = Origami_Elbow_CreasePattern(lengths, ls, n, infostruct(index+3).h1, ...
    infostruct(index+3).h2, theta2, mirror, split, tuckangle);

infostruct(index+3).m = m;
infostruct(index+3).lmax = lmax;
infostruct(index+3).n = n;
infostruct(index+3).type = dataFoldD;
infostruct(index+3).name = "Elbow";
infostruct(index+3).theta = theta2;
infostruct(index+3).dw = r*abs(tan(theta2 / 2));
infostruct(index+3).oc = Op(:, 4) + r*abs(tan(theta1 / 2))*ap + t' + ...
    r*(abs(tan(theta1 / 2))+abs(tan(theta2 / 2)))*tunit';

% If the height of any segment is 0, edit so that lines are not printed
for i = index:index+3
    
    if infostruct(i).lmax == 0
        
        % Replaces x and y values with null, duplicate lines "erased"
        for j = 1:size(infostruct(i).type, 2)
            
            infostruct(i).type(j).x = [];
            infostruct(i).type(j).y = [];
            
        end
        
    end  
    
end

% % Add field for tracking lmaxnet
% infostruct(1).lmaxnet = infostruct(1).lmax;
% 
% msum = 0;
% lmax_sum = 0;
% 
% % figure()
% set(gcf, 'color', 'w')
% hold on
% 
% % Loop through indices to plot
% for i = index:index+3
%     
%     [msum, lmax_sum] = DataFoldAppend(infostruct(i).type, ...
%     infostruct, i, msum, lmax_sum);
%     
% end
% 
% msum = 0;
% lmax_sum = 0;
% 
% % Loop through indices to plot duplication
% for i = index:index+3
%     
%     [msum, lmax_sum] = DataFoldDuplicate(infostruct(i).type, ...
%         infostruct, i, msum, lmax_sum, 'triple');
%     
% end
% 
% for i = 1:size(infostruct, 2)
% 
%     if i == 1
%         lmaxtotal = infostruct(i).lmax;
%     else
%         lmaxtotal = infostruct(i).lmax + lmaxtotal;
%     end
% 
% end
% 
% % Figure adjustments
% daspect([1, 1, 1])
% axis off
% 
% close all


end
