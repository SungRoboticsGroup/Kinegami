function [TransformStruct, JointStructNew, Nnew] = JointPlacementB(D, r, n, JointStruct, N, theta_mod, fingertip, plotoption)
% JOINTPLACEMENTB - This algorithm identifies the joint location
% iteratively such that the new joint location is always 4r apart from all
% the previous joints. In addition, waypoints are introduced to reroute the
% centerline of the link to avoid path intersection. This algorithm shows
% the existence of non-self-intersecting Dubins-specified links.

% Inputs:
%   D               - the D-H parameter table of values: i x [a, alpha, d,
%                     theta].
%   r               - desired radius of folded origami linkage.
%   n               - number of sides of folded origami linkage.
%   JointStruct     - a data structure that contains information about
%                     joint parameters, frames, and connection pathways.
%   N               - variable used for iteration. Equivalent to: (number
%                     of joints of system) - 1.
%   theta_mod       - revolute joint parameter for use within
%                     RotationalMatrix.m.
%   fingertip       - string input ('x', 'y', or 'z') used for fingertip
%                     orientation assignment.
%   plotoption      - string input which dicates plotting.

% Outputs:
%   TransformStruct - updated data structure which contains information
%                     about bounding spheres, associated planes, and other
%                     related information.
%   JointStructNew  - updated data structure which additionally contains
%                     information about waypoint data.
%   Nnew            - value which specifies the size of JointStructNew.
%                     Defined within this function.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Wei-Hsi Chen <weicc@seas.upenn.edu>
% Last Edited 1/25/2023
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% Does this process remain for JointPlacement? Unsure but kept for now.
for i = (N+1):-1:1
    
    % The i index refers to the lower right value for regular and upper
    % left value for inverse
    TransformStruct(i).T = HomogeneousTransform(i, D);
    TransformStruct(i).inverse = InverseHomogeneousTransform(i, D);
    
    % This will be used to find the value of 0{T}N+1
    if i == (N+1)
        TransformStruct(N+1).net = TransformStruct(i).T;
    else
        TransformStruct(N+1).net = TransformStruct(i).T * TransformStruct(N+1).net;
    end
    
end

% Extract O(N+1) and Oc(N+1) from NetTransform (4x4 matrix)
Ox = TransformStruct(N+1).net(1:3, 1);
Oy = TransformStruct(N+1).net(1:3, 2);
Oz = TransformStruct(N+1).net(1:3, 3);
Oc = TransformStruct(N+1).net(1:3, 4);

% O and Oc, initial axis of translation
TransformStruct(N+1).O = [Ox, Oy, Oz, Oc];

if strcmp(fingertip, 'z')
    TransformStruct(N+1).Oc = [Oz, Ox, Oy, Oc];
elseif strcmp(fingertip, 'x')
    TransformStruct(N+1).Oc = [Ox, Oy, Oz, Oc];
elseif strcmp(fingertip, 'y')
    TransformStruct(N+1).Oc = [Oy, Oz, Ox, Oc];
end

TransformStruct(N+1).zglobal = Oz.';
TransformStruct(N+1).xglobal = Ox.';

% Define r and oi for initial sphere (recall .oi etc. is row)
TransformStruct(N+1).oi = Oc.';
TransformStruct(N+1).rs = r;

% Assume value for beta to be constant
beta = pi/4; % [rad]

for i = 1:N+1
    
    % Assign r fields for sphere calculations
    
    % If revolute or extended revolute
    if JointStruct(i).type == 'R'
        
        TransformStruct(i).rs = r*sin(((n - 2)*pi) / (2*n))* ...
            tan(JointStruct(i).qm/ 4);
        
    % If extended revolute
    elseif JointStruct(i).type == 'E'
        
        TransformStruct(i).rs = (r*sin(((n - 2)*pi) / (2*n))* ...
            tan(JointStruct(i).qm/ 4)) + ((JointStruct(i).h1+ ...
            JointStruct(i).h2)/2);    
        
    % If prismatic
    elseif JointStruct(i).type == 'P'
        
        TransformStruct(i).rs = 1/4*JointStruct(i).q0*(2 + csc(beta));
    
    % If waypoint
    elseif JointStruct(i).type == 'W'
        
        TransformStruct(i).rs = 0;
    
    % Otherwise, for fingetip
    else
        
        TransformStruct(i).rs = r*sin(((n - 2)*pi) / (2*n))* ...
            tan(JointStruct(i).qm/4);
        
    end
    
end

% Loop through homogeneous transforms
for i = N:-1:1
    
    % Multiply inverse matrix and standard homogeneous to compute net
    TransformStruct(i).net = TransformStruct(i+1).net * ...
        TransformStruct(i+1).inverse;
    
    % Extract Oi from TransformStruct(i).net
    Ox = TransformStruct(i).net(1:3, 1);
    Oy = TransformStruct(i).net(1:3, 2);
    Oz = TransformStruct(i).net(1:3, 3);
    Oc = TransformStruct(i).net(1:3, 4);
    
    % If revolute or extended revolute
    if JointStruct(i).type == 'R' || JointStruct(i).type == 'E'
        
        Ox = RotationalMatrix(Oz, theta_mod(i))*Ox;
        Oy = RotationalMatrix(Oz, theta_mod(i))*Oy;
        
    end

    % O Computation
    TransformStruct(i).O = [Ox, Oy, Oz, Oc];
    
    % Translation Axis
    TransformStruct(i).zglobal = Oz.';
    
    % X Axis
    TransformStruct(i).xglobal = Ox.';
    
    % Define r and oi for spheres
    TransformStruct(i).oi = Oc.';

end

% Initial visualization of spheres
if strcmp(plotoption, 'on') == 1
    figure()
    set(gcf, 'color', 'w')
    hold on
    grid on
end

% Initial Sphere Visualization (will plot if plotoption is true)
for i = 1:N+1
    
    [TransformStruct(i).oidata] = SphericalSampling(TransformStruct(i).oi, ...
        TransformStruct(i).rs, 'none', plotoption);
    
end

% Initial bounding sphere is same as first sphere
TransformStruct(N+1).oib = TransformStruct(N+1).oi;
TransformStruct(N+1).rb = TransformStruct(N+1).rs;
TransformStruct(N+1).boundingdata = TransformStruct(N+1).oidata;

% Assignment Looping
for i = N+1:-1:2
    
    % Normal vector adjustment (column vector)
    TransformStruct(i).normal = -TransformStruct(i).Oc(:, 1);

    % Plane generation for tangent plane
    TransformStruct(i).tangent(:, 4) = TransformStruct(i).oib.' + ...
      TransformStruct(i).rb * TransformStruct(i).normal;

    TransformStruct(i).tangent(:, 1) = TransformStruct(i).normal;

    % Storing plane data
    TransformStruct(i).P1 = zeros(3, 2);
    TransformStruct(i).P1(:, 1) = TransformStruct(i).tangent(:, 1);
    TransformStruct(i).P1(:, 2) = TransformStruct(i).tangent(:, 4);  

    % Plane generation for parallel plane
    TransformStruct(i).parallel(:, 4) = TransformStruct(i).oib.' + ...
      (TransformStruct(i).rb + TransformStruct(i-1).rs + 4*r) * ...
      TransformStruct(i).normal;

    TransformStruct(i).parallel(:, 1) = TransformStruct(i).normal;

    % Storing parallel plane data
    TransformStruct(i).P1par = zeros(3, 2);
    TransformStruct(i).P1par(:, 1) = TransformStruct(i).parallel(:, 1);
    TransformStruct(i).P1par(:, 2) = TransformStruct(i).parallel(:, 4);  

    % Determine location of first waypoint
    [waypoint1] = IntersectionSolver(TransformStruct(i).P1, ...
      TransformStruct(i).Oc(:, 4).', TransformStruct(i).normal);

    % Store to structure as waypoint1
    TransformStruct(i).waypoint1 = [TransformStruct(i).Oc(:, 1), ...
      TransformStruct(i).Oc(:, 2), TransformStruct(i).Oc(:, 3), ...
      waypoint1];

    % Determine intersection of preceding z-axis (for sphere being placed)
    % and plane, if it occurs.

    % Run PlaneCheck, which will determine if the z-axis is on the far or
    % near side of the parallel plane.
      planeval = PlaneCheck(TransformStruct(i).P1par, TransformStruct(i-1).oi);   

    % For parallel, conflicting case
    % There could be numeric error, so check for "statment < 0.00001"
    % instead of "statement == 0".
    if abs(dot(TransformStruct(i).normal, TransformStruct(i-1).zglobal)) < 0.00001 && ...
            planeval < 0

        % Assign normal vector
        TransformStruct(i).normal2 = -TransformStruct(i-1).zglobal.';

        % Assign new planes
        TransformStruct(i).tangent2(:, 4) = TransformStruct(i).oib.' + ...
            TransformStruct(i).rb * TransformStruct(i).normal2;

        TransformStruct(i).tangent2(:, 1) = TransformStruct(i).normal2;

        % Storing plane data
        TransformStruct(i).P2 = zeros(3, 2);
        TransformStruct(i).P2(:, 1) = TransformStruct(i).tangent2(:, 1);
        TransformStruct(i).P2(:, 2) = TransformStruct(i).tangent2(:, 4);  

        % Plane generation for parallel plane
        TransformStruct(i).parallel2(:, 4) = TransformStruct(i).oib.' + ...
            (TransformStruct(i).rb + TransformStruct(i-1).rs + 4*r) * ...
            TransformStruct(i).normal2;

        TransformStruct(i).parallel2(:, 1) = TransformStruct(i).normal2;

        % Storing parallel plane data in useful format
        TransformStruct(i).P2par = zeros(3, 2);
        TransformStruct(i).P2par(:, 1) = TransformStruct(i).parallel2(:, 1);
        TransformStruct(i).P2par(:, 2) = TransformStruct(i).parallel2(:, 4);  

        % Determining waypoint2
        ref2 = waypoint1 + r*TransformStruct(i).normal; 
        [waypoint2] = IntersectionSolver(TransformStruct(i).P2, ...
            ref2, TransformStruct(i).normal2);

        % Rotational Frame Transform
        rotationvec = cross(TransformStruct(i).normal, TransformStruct(i).normal2);
        rotationvec = rotationvec/norm(rotationvec);
      
        [rotation_matrix] = RotationalMatrix(rotationvec, pi/2);
      
        OcMatrix = [TransformStruct(i).Oc(:, 1), TransformStruct(i).Oc(:, 2), ...
            TransformStruct(i).Oc(:, 3)];
        ORotation =  rotation_matrix * OcMatrix;
      
        % Store to structure as waypoint2
        TransformStruct(i).waypoint2 = [ORotation(:, 1), ORotation(:, 2), ...
            ORotation(:, 3), waypoint2];

    % For non-parallel case, non-conflicting case
    else
        
        % Assign n2 vector
        TransformStruct(i).normal2 = TransformStruct(i).normal;

        % Assign new planes
        TransformStruct(i).P2 = TransformStruct(i).P1;
        TransformStruct(i).P2par = TransformStruct(i).P1par;

        % Assign waypoint2
        TransformStruct(i).waypoint2 = TransformStruct(i).waypoint1;

    end

    % Sphere bounding process

    % Oc assignment
    % Joint Type Dictates Ox, Oy, Oz, Oc order. Clearly define these.
    Ox = TransformStruct(i-1).O(:, 1);
    Oy = TransformStruct(i-1).O(:, 2);
    Oz = TransformStruct(i-1).O(:, 3);
    Oc = TransformStruct(i-1).O(:, 4);
    
    % Since fingertip case (N+1) is already addressed in lines 29-44,
    % worry only about Prismatic and Revolute cases
    if JointStruct(i-1).type == 'R' || JointStruct(i).type == 'E'
        
        TransformStruct(i-1).Oc = [Ox, Oy, Oz, Oc];
        
    elseif JointStruct(i-1).type == 'P'
        
        TransformStruct(i-1).Oc = [Oz, Ox, Oy, Oc];        
                
    end
    
    % Now that we have bare bones template for Oc, vary values
    % Begin by determining u and t    
    if -TransformStruct(i-1).Oc(:, 1).' * TransformStruct(i).normal2 >= 0
        
        u = 1;       
        
    elseif -TransformStruct(i-1).Oc(:, 1).' * TransformStruct(i).normal2 < 0
        
        u = -1;
        
    end
    
    % Optimization based on constraints to determine value of t

    % Statement for minimization
    obj = @(delta) norm(Oc + delta*Oz - ...
        TransformStruct(i).waypoint2(:, 4));

    % Initial Guess
    delta0 = 1;
    
    % Run PlaneFind to determine variables for a, b, c, d for plane in
    % question.    
    [a, b, c, d] = PlaneFind(TransformStruct(i).P2par);

    % Inequality requirement:
    % a*oi(1) + b*oi(2) + c*oi(3) + d >= 0

    % Express other parameters as empty cells (not used)
    A = -dot([a b c], Oz');
    B = dot([a b c], Oc')+d;

    % To avoid potential numeric error:
    if abs(A) < 0.00001

        A = 0;
        
    end    
    
    % Output optimal value for delta (or t)
    t_val = fmincon(obj, delta0, A, B);
    
    % Now, use values of u and t to find Oc(i-1)
    if JointStruct(i-1).type == 'R' || JointStruct(i-1).type == 'E'
        
        TransformStruct(i-1).Oc = [u*Ox, u*Oy, Oz, Oc + t_val*Oz];
        
    elseif JointStruct(i-1).type == 'P'
        
        TransformStruct(i-1).Oc = [u*Oz, u*Ox, Oy, Oc + t_val*Oz];        
                
    end
    
    
    % New sphere creation
    [TransformStruct(i-1).oidata] = SphericalSampling(TransformStruct(i-1).Oc(:, 4), ...
        TransformStruct(i-1).rs, 'none', plotoption);
    
    % Bounding sphere data
    TransformStruct(i-1).boundingdata = [TransformStruct(i).boundingdata; ...
        TransformStruct(i-1).oidata];
    
    % Find new oib and rb values
    [R, C, Xb] = ExactMinBoundSphere3D(TransformStruct(i-1).boundingdata);
    
    TransformStruct(i-1).oib = C;
    TransformStruct(i-1).rb = R;
    TransformStruct(i-1).Xb = Xb;
    
    % Output plot of new bounding sphere
    [TransformStruct(i-1).boundingsummary] = SphericalSampling(TransformStruct(i-1).oib, ...
        TransformStruct(i-1).rb, 'none', plotoption);
    
end

% Create O DataSet (Line 21 of Script)
% Create new JointStruct data structure which includes waypoints
Nnew = 3*(N);
JointStructNew(Nnew+1) = struct();

% JointStructNew(1) is the same as JointStruct(1), then JointStruct(2)
% corresponds to JointStructNew(4), etc.
index = 0;
for i = 1:3:(Nnew+1)
    
    index = index+1;
    JointStructNew(i).qm = JointStruct(index).qm;
    JointStructNew(i).q0 = JointStruct(index).q0;
    JointStructNew(i).type = JointStruct(index).type;
    JointStructNew(i).nz = JointStruct(index).nz;
    
end

% Populate waypoint values
index = 1;
for i = 2:3:(Nnew)
    
    index = index + 1;
    % Waypoint 1
    JointStructNew(i).Op = TransformStruct(index).waypoint2;
    JointStructNew(i).Od = TransformStruct(index).waypoint2;
    JointStructNew(i).type = 'W';
    
    % Waypoint 2
    JointStructNew(i+1).Op = TransformStruct(index).waypoint1;
    JointStructNew(i+1).Od = TransformStruct(index).waypoint1;
    JointStructNew(i+1).type = 'W';
    
end

% Final sphere ("joint") visualization
if strcmp(plotoption, 'on') == 1
    figure()
    set(gcf, 'color', 'w')
    hold on
end

black = [0,0,0];
red = [0.6350, 0.0780, 0.1840];
orange = [0.8500, 0.3250, 0.0980];
yellow = [0.929, 0.694, 0.125];
green = [0.466, 0.674, 0.188];
blue = [0, 0.447, 0.741];
purple = [0.4940, 0.1840, 0.5560];

colorvector = [black; red; orange; yellow; green; blue; purple];

for i = 1:N+1
        
    % Plot final new sphere locations
    [TransformStruct(i).demonew] = SphericalSampling(TransformStruct(i).Oc(:, 4), ...
        TransformStruct(i).rs, colorvector(i, :), plotoption); 
       
    % Plot vectors which demonstrate the axis along which each sphere is
    % restricted to move
    if strcmp(plotoption, 'on') == 1
        quiver3(TransformStruct(i).Oc(1, 4), TransformStruct(i).Oc(2, 4), ...
            TransformStruct(i).Oc(3, 4), TransformStruct(i).zglobal(1), ...
            TransformStruct(i).zglobal(2), TransformStruct(i).zglobal(3), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 0.8*colorvector(i, :));

        quiver3(TransformStruct(i).Oc(1, 4), TransformStruct(i).Oc(2, 4), ...
            TransformStruct(i).Oc(3, 4), TransformStruct(i).xglobal(1), ...
            TransformStruct(i).xglobal(2), TransformStruct(i).xglobal(3), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 'k');

        grid on
    end
end

% Plot lines connecting consecutive spheres

% Initialize
xcenters = zeros(N+1, 1);
ycenters = zeros(N+1, 1);
zcenters = zeros(N+1, 1);

% Loop through to create centerpoint vectors
for i = 1:N+1
    
    xcenters(i) = TransformStruct(i).Oc(1, 4);
    ycenters(i) = TransformStruct(i).Oc(2, 4);
    zcenters(i) = TransformStruct(i).Oc(3, 4);   
    
end

% Plot
if strcmp(plotoption, 'on') == 1
    plot3(xcenters(:, 1), ycenters(:, 1), zcenters(:, 1), 'Color', 'k', ...
        'Linewidth', 4)
end

% We can now use the list of TransformStruct.Oc to provide us information
% on the new centroids of all of these spheres (joints)

% New visualization
if strcmp(plotoption, 'on') == 1
    hold off
    figure()
    set(gcf, 'color', 'w')
    hold on
end

for i = 1:N+1  
    
    [TransformStruct(i).demonew] = SphericalSampling(TransformStruct(i).Oc(:, 4), ...
    TransformStruct(i).rs, colorvector(i, :), plotoption); 
    
    if strcmp(plotoption, 'on') == 1
        % Plot vectors which demonstrate the new Oc 
        quiver3(TransformStruct(i).Oc(1,4), TransformStruct(i).Oc(2,4), ...
            TransformStruct(i).Oc(3,4), TransformStruct(i).Oc(1,1), ...
            TransformStruct(i).Oc(2,1), TransformStruct(i).Oc(3,1), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 'red');

        quiver3(TransformStruct(i).Oc(1,4), TransformStruct(i).Oc(2,4), ...
            TransformStruct(i).Oc(3,4), TransformStruct(i).Oc(1,2), ...
            TransformStruct(i).Oc(2,2), TransformStruct(i).Oc(3,2), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 'green');

        quiver3(TransformStruct(i).Oc(1,4), TransformStruct(i).Oc(2,4), ...
            TransformStruct(i).Oc(3,4), TransformStruct(i).Oc(1,3), ...
            TransformStruct(i).Oc(2,3), TransformStruct(i).Oc(3,3), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 'blue');

        grid on
    end
    
end

% Plot
if strcmp(plotoption, 'on') == 1
    plot3(xcenters(:, 1), ycenters(:, 1), zcenters(:, 1), 'Color', 'k', ...
        'Linewidth', 4)
end

end