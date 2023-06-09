% KINEGAMI_TEMPLATE - Generates a crease pattern for the given robot
% specification. Fill in the Kinematic parameters '??' and specify the user
% options. Written in the "selfassign" option.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Wei-Hsi Chen <weicc@seas.upenn.edu>
% Daniel Feshbach <feshbach@seas.upenn.edu>
% Last Edited 04/25/2023
%
% Copyright (C) 2023 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


clear
close all
clc

% Add all the folders and subfolders to the search path
addpath(genpath(fileparts(mfilename('fullpath'))));

%% USER OPTIONS - Change Prior to Running (if necessary)

% Determines which joint placement method will be implemented:
% For the general method, use JointPlacementA.m ('placementA'); 
% To guarantee no self-intersection, use JointPlacementB.m ('placementB');
% To locate the joints manually, use SelfAssign.m ('selfassign')
jointselect = 'selfassign';


% Determines whether the user wishes their elbow fittings to have visible
% tucks ('on' - recommended) or appear with only the lower outlines ('off')
elbow_tuck = 'on';

% Specify whether elbow splitting should occur past pi/2 ('on'/'off')
split = 'on';

% Determines whether the user wishes to print 3 iterations of the print
% pattern ('triple' - recommended) or 2 ('double')
triple = 'triple';

% Specify whether DXF generation and save file should occur ('on'/'off')
DXF = 'on';

% Specify whether or not pre-segmentation (for printing) is enabled
% ('on'/'off')
segmentation = 'off';

% Specify whether intermediary plots should be run ('on'/'off'). 'off' is
% recommended for faster computational time, 'on' is recommended for more
% in-depth analysis.
plotoption = 'off';

% Specify whether initial tube plotting (prior to first joint location)
% should be included. ('on'/'off')
tubeinit = 'off';

%% KINEMATIC CHAIN SPECIFICATION

% Specify number of sides and the circumradius [m] for the polygon base 
% of the prism tube
nsides = 4;
r = 0.04;

% Input the kinematic chain robot specifications
% Number of joints and initialize DH Parameter table D
% n equals to the number of joints plus one fingertip
n = 8;
D = zeros(n,4);

% Specify DH Parameters nX4, in the order of "Link length (a)", 
% "Link twist (alpha)", "Joint offset (d)", and "Joint angle (theta)".
% D = [??, ??, ??, ??;...
%     ...
%      ??, ??, ??, ??];

% Specify joint information as a row vector:
% Contains n elements for each vector
% Types of joints: 'R': Revolute joint, 'E': Extended Revolute Joint, ...
% 'P': Prismatic joint, 'F': Fingertip, 'W': Waypoint (not a joint)
TYPE = ['W', 'E', 'E', 'W', 'W', 'E', 'E', 'W']; 

% Maximum joint range (row vec.)
Qm = [pi/2, pi, pi, pi/2, pi/2, pi, pi, pi/2]; 

% Initial joint configuration (row vec.)
Q0 = [0, 0, 0, 0, 0, 0, 0, 0]; 

% Specify the angle modification utilized (recommended: zeros(n)) (row vec.)
theta_mod = [0, 0, 0, 0, 0, 0, 0, 0]; 

% Specify the z-axis modifications
% To leave the z-axis placement up to the algorithm, use NaN
% z_mod = [0   0   0   0   0   0   0   0];
z_mod = [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];
% Solution found when z_mod is all NaN:



% Layer of recursive sink gadget for revolute joint (row vec.)
Nz = [1, 1, 1, 1, 1, 1, 1, 1]; 

% Specify the orientation of the fingertip: 'x', 'y', or 'z'
fingertip = 'x';

% Initialize JointStruct
JointStruct(n) = struct();
for i = 1:n
    JointStruct(i).qm = Qm(i);
    JointStruct(i).q0 = Q0(i);
    JointStruct(i).type = TYPE(i);
    JointStruct(i).nz = Nz(i);
end
N = size(JointStruct, 2) - 1;
TransformStruct(N+1) = struct();

%% SELF-ASSIGNED JOINT LOCATION

% If the selfassign tag is applied, provide Oc for each joint
if strcmp(jointselect, 'selfassign') == 1
    
    % Specify the orientation of the fingertip: 'x', 'y', or 'z'
    fingertip = 'x';
    
    % Lengths [m]
    l1 = 0.2;
    l2 = 0.2;
    d = 2*l1;
    
    % Specify the frame (3X4) of each individul joint in global frame
    TransformStruct(1).Oc = [1, 0, 0, -l2; ...
                             0, 1, 0, 0; ...
                             0, 0, 1, 0];
                         
    TransformStruct(2).Oc = [1, 0, 0, 0; ...
                             0, 1, 0, 0; ...
                             0, 0, 1, 0];
                         
    TransformStruct(3).Oc = [cos(-pi/3), -sin(-pi/3), 0, l2*cos(-pi/3); ...
                             sin(-pi/3), cos(-pi/3), 0, l2*sin(-pi/3); ...
                             0, 0, 1, 0];  
                     
    TransformStruct(4).Oc = [cos(-2*pi/3), -sin(-2*pi/3), 0, l2*cos(-pi/3) - ((d*cos(pi/3))/3); ...
                             sin(-2*pi/3), cos(-2*pi/3), 0, l2*sin(-pi/3) - ((d*sin(pi/3)/3)); ...
                             0, 0, 1, -0.125];   
                         
    TransformStruct(5).Oc = [cos(-2*pi/3), -sin(-2*pi/3), 0, l2*cos(-pi/3) - (2*(d*cos(pi/3))/3); ...
                             sin(-2*pi/3), cos(-2*pi/3), 0, l2*sin(-pi/3) - (2*(d*sin(pi/3))/3); ...
                             0, 0, 1, -0.125];   
                         
    TransformStruct(6).Oc = [cos(2*pi/3), -sin(2*pi/3), 0, (l2*cos(-pi/3)) - d*cos(pi/3); ...
                             sin(2*pi/3), cos(2*pi/3), 0, (l2*sin(-pi/3)) - d*sin(pi/3); ...
                             0, 0, 1, 0];                        
                                               
    TransformStruct(7).Oc = [-1, 0, 0, (l2*cos(-pi/3)) - d*cos(pi/3) + (l2*cos(2*pi/3)); ...
                             0, -1, 0, (l2*sin(-pi/3)) - d*sin(pi/3) + (l2*sin(2*pi/3)); ...
                             0, 0, 1, 0];  

    TransformStruct(8).Oc = [-1, 0, 0, (l2*cos(-pi/3)) - d*cos(pi/3) + (l2*cos(2*pi/3)) - l1; ...
                             0, -1, 0, (l2*sin(-pi/3)) - d*sin(pi/3) + (l2*sin(2*pi/3)); ...
                             0, 0, 1, 0];                           
                         
else  
    % Otherwise, do nothing new
end

%% RUN KINEGAMI

% Run Kinegami code
[infostruct, TransformStruct, DataNet, JointStruct] = Kinegami(D, r, nsides, JointStruct, ...
    elbow_tuck, triple, theta_mod, z_mod, fingertip, TransformStruct, ...
    DXF, split, segmentation, plotoption, jointselect, tubeinit);
