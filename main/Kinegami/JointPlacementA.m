function [TransformStruct] = JointPlacementA(D, r, n, JointStruct, N, theta_mod, z_mod, fingertip, plotoption)
% JOINTPLACEMENTA - This algorithm identifies the joint location
% iteratively such that the new joint location is always 4r apart from all
% the previous joints. This algorithm shows the existence of
% Dubins-specified links.

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

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Wei-Hsi Chen <weicc@seas.upenn.edu>
% Daniel Feshbach <feshbach@seas.upenn.edu>
% Last Edited 5/11/2023
%
% Copyright (C) 2023 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.

% Assume value for beta to be constant
beta = pi/4; % [rad]

for i = 1:(N+1)
    
    if isfield(JointStruct, 'zAxisRelative')
        TransformStruct(i).T = HomogeneousTransformOZ( ...
            JointStruct(i).oRelative, JointStruct(i).zAxisRelative);
    else
        TransformStruct(i).T = HomogeneousTransformDH_KD(i, D);
    end
    
    if i == 1
        TransformStruct(1).net = TransformStruct(1).T;
    else
        TransformStruct(i).net = TransformStruct(i-1).net * TransformStruct(i).T;
    end
    
    % Extract Oi from TransformStruct(i).net
    Ox = TransformStruct(i).net(1:3, 1);
    Oy = TransformStruct(i).net(1:3, 2);
    Oz = TransformStruct(i).net(1:3, 3);
    Oc = TransformStruct(i).net(1:3, 4);

    % O Computation
    TransformStruct(i).O = [Ox, Oy, Oz, Oc];
    
    if i == N+1
        if strcmp(fingertip, 'z')
            TransformStruct(N+1).Oc = [Oz, Ox, Oy, Oc];
        elseif strcmp(fingertip, 'x')
            TransformStruct(N+1).Oc = [Ox, Oy, Oz, Oc];
        elseif strcmp(fingertip, 'y')
            TransformStruct(N+1).Oc = [Oy, Oz, Ox, Oc];
        end
        %TransformStruct(N+1).rs = r; % I think this ends up getting over-written
    end

    % Translation Axis
    TransformStruct(i).zHatAbsolute = Oz.';
    
    % X Axis
    TransformStruct(i).xHatAbsolute = Ox.';
    
    % Define r and oi for spheres
    TransformStruct(i).oi = Oc.';


    % Assign r fields for sphere calculations
    
    % If revolute
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
    
    % Otherwise, for fingetip
    else
        
        TransformStruct(i).rs = r*sin(((n - 2)*pi) / (2*n))* ...
            tan(JointStruct(i).qm/4);
        
    end
    
end

black = [0,0,0];
red = [0.6350, 0.0780, 0.1840];
orange = [0.8500, 0.3250, 0.0980];
yellow = [0.929, 0.694, 0.125];
green = [0.466, 0.674, 0.188];
blue = [0, 0.447, 0.741];
purple = [0.4940, 0.1840, 0.5560];

colorvector = [black; red; orange; yellow; green; blue; purple];

% Initial visualization of spheres
if strcmp(plotoption, 'on') == 1
    figure()
    set(gcf, 'color', 'w')
    hold on
end

for i = 1:N+1
    
    [TransformStruct(i).demo] = SphericalSampling(TransformStruct(i).oi, ...
        TransformStruct(i).rs, 'none', plotoption); 

    % Plot vectors which demonstrate the axis along which each sphere is
    % restricted to move
    if strcmp(plotoption, 'on') == 1
        z_scale = 10;

        quiver3(TransformStruct(i).oi(1), TransformStruct(i).oi(2), ...
            TransformStruct(i).oi(3), z_scale*TransformStruct(i).zHatAbsolute(1), ...
            z_scale*TransformStruct(i).zHatAbsolute(2), z_scale*TransformStruct(i).zHatAbsolute(3), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 0.8*colorvector(i, :));

        % quiver3(TransformStruct(i).oi(1), TransformStruct(i).oi(2), ...
        %     TransformStruct(i).oi(3), -z_scale*TransformStruct(i).zHatAbsolute(1), ...
        %     -z_scale*TransformStruct(i).zHatAbsolute(2), -z_scale*TransformStruct(i).zHatAbsolute(3), ...
        %     'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 0.8*colorvector(i, :));

        quiver3(TransformStruct(i).oi(1), TransformStruct(i).oi(2), ...
            TransformStruct(i).oi(3), TransformStruct(i).xHatAbsolute(1), ...
            TransformStruct(i).xHatAbsolute(2), TransformStruct(i).xHatAbsolute(3), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 'k');

        grid on
    end
end

% Sphere bounding fundamentals for (N+1) ---------------------------

% Sphere bounding data
if strcmp(plotoption, 'on') == 1
    figure()
    set(gcf, 'color', 'w')
    hold on
end

for i=1:N+1
    TransformStruct(i).delta_z = z_mod(i);
end

if isnan(TransformStruct(N+1).delta_z)
    TransformStruct(N+1).delta_z = 0;
end
TransformStruct(N+1).oinew = TransformStruct(N+1).oi + ...
    TransformStruct(N+1).delta_z * TransformStruct(N+1).zHatAbsolute;

[TransformStruct(N+1).demo] = SphericalSampling(TransformStruct(N+1).oinew, ...
    TransformStruct(N+1).rs, 'none', plotoption);

% Concatenate to express array of spheres
TransformStruct(N+1).democumul = TransformStruct(N+1).demo;

% Find new sphere center and store to index-1.optimized
[R, C, Xb] = ExactMinBoundSphere3D(TransformStruct(N+1).democumul);

TransformStruct(N+1).Xb = Xb;
TransformStruct(N+1).oinew = C;
TransformStruct(N+1).rnew = R;

% Output new plot of cumulative sphere
[TransformStruct(N+1).democumul] = SphericalSampling(TransformStruct(N+1).oinew, ...
    TransformStruct(N+1).rnew, 'none', plotoption);

TransformStruct(N+1).oilarge = TransformStruct(N+1).oi;

% Loop through sphere bounding and minimization (Fun!)
for i = (N+1):-1:2 
   [TransformStruct] = SphereMinimization(TransformStruct, i, r, 'none', plotoption);  
end

% Final sphere ("joint") visualization
if strcmp(plotoption, 'on') == 1
    figure()
    set(gcf, 'color', 'w')
    hold on
end

for i = 1:N+1
    
    % Original spheres
%     [TransformStruct(i).demo] = SphericalSampling(TransformStruct(i).oi, ...
%     TransformStruct(i).r, 0.7*colorvector(i, :)); 
% 
%     quiver3(TransformStruct(i).oi(1), TransformStruct(i).oi(2), ...
%         TransformStruct(i).oi(3), TransformStruct(i).axis(1), ...
%         TransformStruct(i).axis(2), TransformStruct(i).axis(3), ...
%         'AutoScaleFactor', 1, 'Linewidth', 1)
    
    % Plot final new sphere locations
    i
    [TransformStruct(i).demonew] = SphericalSampling(TransformStruct(i).oinew, ...
        TransformStruct(i).rs, colorvector(i, :), plotoption); 
    
    
    % Plot vectors which demonstrate the axis along which each sphere is
    % restricted to move
    if strcmp(plotoption, 'on') == 1
        
        % The z axes of the original joint origins, for reference
        quiver3(TransformStruct(i).oi(1), TransformStruct(i).oi(2), ...
            TransformStruct(i).oi(3), TransformStruct(i).zHatAbsolute(1), ...
            TransformStruct(i).zHatAbsolute(2), TransformStruct(i).zHatAbsolute(3), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 3, 'Color', 0.8*colorvector(i, :));        
        
        quiver3(TransformStruct(i).oinew(1), TransformStruct(i).oinew(2), ...
            TransformStruct(i).oinew(3), TransformStruct(i).zHatAbsolute(1), ...
            TransformStruct(i).zHatAbsolute(2), TransformStruct(i).zHatAbsolute(3), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 0.8*colorvector(i, :));

        quiver3(TransformStruct(i).oinew(1), TransformStruct(i).oinew(2), ...
            TransformStruct(i).oinew(3), TransformStruct(i).xHatAbsolute(1), ...
            TransformStruct(i).xHatAbsolute(2), TransformStruct(i).xHatAbsolute(3), ...
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
    
    xcenters(i) = TransformStruct(i).oinew(1);
    ycenters(i) = TransformStruct(i).oinew(2);
    zcenters(i) = TransformStruct(i).oinew(3);   
    
end

% Plot
if strcmp(plotoption, 'on') == 1
    plot3(xcenters(:, 1), ycenters(:, 1), zcenters(:, 1), 'Color', 'k', ...
        'Linewidth', 4)
end
    

% Loop through joint coordinate reassignment
for i = 1:N+1
    
    % Redefine x y and z, as well as centroid, based on sphere reassignment
    Ox = TransformStruct(i).net(1:3, 1);
    Oy = TransformStruct(i).net(1:3, 2);
    Oz = TransformStruct(i).net(1:3, 3);
    Oc = TransformStruct(i).oinew.';
    
    % If revolute or extended revolute
    if JointStruct(i).type == 'R' || JointStruct(i).type == 'E'
        
        ai = RotationalMatrix(Oz, theta_mod(i))*Ox;
        bi = RotationalMatrix(Oz, theta_mod(i))*Oy;
        
        TransformStruct(i).Oc = [ai, bi, Oz, Oc];
    
    % If prismatic
    elseif JointStruct(i).type == 'P'
        
        TransformStruct(i).Oc = [Oz, Ox, Oy, Oc];      
    
    % Other joints
    else
        
        % If fingertip assignment is x
        if strcmp(fingertip, 'x') == 1
            
            TransformStruct(i).Oc = [Ox, Oy, Oz, Oc];
        
        % If fingertip assignment is y
        elseif strcmp(fingertip, 'y') == 1
            
            TransformStruct(i).Oc = [Oy, Oz, Ox, Oc];
        
        % If fingertip assignment is z
        else
            
            TransformStruct(i).Oc = [Oz, Ox, Oy, Oc];
            
        end
        
    end
    
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
    
    [TransformStruct(i).demonew] = SphericalSampling(TransformStruct(i).oinew, ...
    TransformStruct(i).rs, colorvector(i, :), plotoption); 
    
    if strcmp(plotoption, 'on') == 1
        % Plot vectors which demonstrate the new Oc 
        quiver3(TransformStruct(i).oinew(1), TransformStruct(i).oinew(2), ...
            TransformStruct(i).oinew(3), TransformStruct(i).Oc(1,1), ...
            TransformStruct(i).Oc(2,1), TransformStruct(i).Oc(3,1), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 'red');

        quiver3(TransformStruct(i).oinew(1), TransformStruct(i).oinew(2), ...
            TransformStruct(i).oinew(3), TransformStruct(i).Oc(1,2), ...
            TransformStruct(i).Oc(2,2), TransformStruct(i).Oc(3,2), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 'green');

        quiver3(TransformStruct(i).oinew(1), TransformStruct(i).oinew(2), ...
            TransformStruct(i).oinew(3), TransformStruct(i).Oc(1,3), ...
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