function [tmin, theta1min, theta2min] = dubinsCSC_HG10(r, Od, Op)
% SOLVEDUBINS3D - Find the CSC Dubins path from frame Op to Od.
% Outputs row vector and two radian values. Ensure that Od and Op are
% normalized prior to entry. Only guaranteed to work if op and od are
% at least 4r apart. Implementation of S. Hota and D. Ghose, 
% "Optimal Geometrical Path in 3D with Curvature Constraint", IROS 2010.

% Inputs:
%   r           - desired radius of folded origami linkage.
%   Od          - distal, or terminal, frame of dubins pathing.
%   Op          - proximal, or initial, frame of dubins pathing.

% Outputs:
%   tmin        - row vector associated with the shortest path for Dubins
%                 path linkage.
%   theta1min   - theta value (in radians) for use in RotationalMatrix.m,
%                 for utility in determining rotational matrix around wp.
%   theta2min   - theta value (in radians) for use in RotationalMatrix.m,
%                 for utility in determining rotational matrix around wd.

% Authors: 
% Cynthia Sung <crsung@seas.upenn.edu>
% Lucien Peach <peach@seas.upenn.edu>
% Wei-Hsi Chen <weicc@seas.upenn.edu>
% Lasted edited 7/27/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% input params
od = Od(:, 4).';
ad = Od(:, 1).';
bd = Od(:, 2).';

op = Op(:, 4).';
ap = Op(:, 1).';
bp = Op(:, 2).';

% normalize 'velocity vectors'
ad = ad / norm(ad);
ap = ap / norm(ap);

% normalize rotational reference vectors
bd = bd / norm(bd);
bp = bp / norm(bp);

% solve for optimal t = [tx, ty, tz]
T0 = od-op+randn(size(od));

% Initialize vectors for error and solutions
Tsol = zeros(4, 3);
Terror = zeros(4, 3);

% provides individual solutions 
Tsol(1, :) = fsolve(@dubinspath1, T0); % + ... +
Tsol(2, :) = fsolve(@dubinspath2, T0); % - ... -
Tsol(3, :) = fsolve(@dubinspath3, T0); % + ... -
Tsol(4, :) = fsolve(@dubinspath4, T0); % - ... +

thetamatrix = zeros(4, 2);
ds = zeros(4, 1);
T_hat = zeros(4, 3);

for i = 1:4
    
    T_normalized = norm(Tsol(i, :));
    T_hat(i, :) = Tsol(i, :)/T_normalized;
    
% Find angle using arccos, the range of the output angle is (0, pi)   
%     thetamatrix(i, 1) = acos(dot(ap, T_hat(i, :)));
%     thetamatrix(i, 2) = acos(dot(ad, T_hat(i, :)));

% Find angle using atan2, the range of the output angle is (-pi, pi)   
%     thetamatrix(i, 1) = atan2(norm(cross(ap, T_hat(i, :))),dot(ap, T_hat(i, :)));
%     thetamatrix(i, 2) = atan2(norm(cross(T_hat(i, :), ad)),dot(T_hat(i, :), ad));
    
    thetamatrix(i, 1) = SignedAngle(ap, T_hat(i, :), cross(ap, T_hat(i, :)));
%     thetamatrix(i, 2) = SignedAngle(T_hat(i, :), ad, cross(T_hat(i, :), ad));
    
    % find the point om after a C + S
    wp = cross(ap, T_hat(i, :));
    wp = wp/norm(wp);
    o2c1 = -cross(ap, wp)';
    o2c1 = o2c1/norm(o2c1);
    om = op' + r*(eye(3) - RotationalMatrix(wp, thetamatrix(i, 1)))* o2c1 + Tsol(i, :)';
    % find the vector from om to od, create a normal vector as the axis of
    % rotation
    om2od = od' - om;
    ntemp = cross(T_hat(i, :), om2od);

    thetamatrix(i, 2) = SignedAngle(T_hat(i, :), ad, ntemp);
    
    
    ds(i) = r*abs(thetamatrix(i, 1)) + T_normalized + r*abs(thetamatrix(i, 2));
    
end

[minds, indexds] = min(ds);

disp(minds);

% Output minimum Tunittor as well as minimum theta1 and theta2 values
tmin = Tsol(indexds, :);
theta1min = thetamatrix(indexds, 1);
theta2min = thetamatrix(indexds, 2);

% double check the error is low
Terror(1, :) = dubinspath1(Tsol(1, :));
Terror(2, :) = dubinspath2(Tsol(2, :));
Terror(3, :) = dubinspath3(Tsol(3, :));
Terror(4, :) = dubinspath4(Tsol(4, :));

% Display solutions and error
disp(Tsol)
disp(Terror)

% Plotting
% --------------------------------------------------------------------

% decompose path in relevant points for plotting
[Cd1, Pd1, Cp1, Pp1] = decomposePath1(Tsol(1,:));

% plotting
hold on
subplot(2, 2, 1)
plotCirclePath(op, od, Cp1, Pp1, Cd1, Pd1, ap, ad, 'b', 'r');

axis equal

% decompose path in relevant points for plotting
[Cd2, Pd2, Cp2, Pp2] = decomposePath2(Tsol(2,:));

% plotting
hold on
subplot(2, 2, 2)
plotCirclePath(op, od, Cp2, Pp2, Cd2, Pd2, ap, ad, 'b', 'r')

axis equal

% decompose path in relevant points for plotting
[Cd3, Pd3, Cp3, Pp3] = decomposePath3(Tsol(3,:));

% plotting
hold on
subplot(2, 2, 3)
plotCirclePath(op, od, Cp3, Pp3, Cd3, Pd3, ap, ad, 'b', 'r')

axis equal

% decompose path in relevant points for plotting
[Cd4, Pd4, Cp4, Pp4] = decomposePath4(Tsol(4,:));

% plotting
hold on
subplot(2, 2, 4)
plotCirclePath(op, od, Cp4, Pp4, Cd4, Pd4, ap, ad, 'b', 'r')

axis equal

AllPathSpecs = zeros(4,4);
[th11,g11,th21,g21] = pathSpecs(Cd1, Pd1, Cp1, Pp1);
AllPathSpecs(1,:) = [th11,g11,th21,g21];
[th12,g12,th22,g22] = pathSpecs(Cd2, Pd2, Cp2, Pp2);
AllPathSpecs(2,:) = [th12,g12,th22,g22];
[th13,g13,th23,g23] = pathSpecs(Cd3, Pd3, Cp3, Pp3);
AllPathSpecs(3,:) = [th13,g13,th23,g23];
[th14,g14,th24,g24] = pathSpecs(Cd4, Pd4, Cp4, Pp4);
AllPathSpecs(4,:) = [th14,g14,th24,g24];
AllPathSpecs

% Functions
% ---------------------------------------------------------------------

    function [Cd, Pd, Cp, Pp] = decomposePath1(T) % + ... +
        % find relevant points based on T vector
        Tnorm = norm(T);
        Tunit = T/Tnorm;
        
        wp = cross(ap, cross(Tunit, ap));
        yp = cross(Tunit, cross(Tunit, ap));
        Cp = op + r * wp / norm(wp); % center of starting circle
        Pp = Cp - r * yp / norm(yp); % point of leaving starting circle
        
        wd = -cross(ad, cross(Tunit, ad));
        yd = -cross(Tunit, cross(Tunit, ad));
        Cd = od + r * wd / norm(wd); % center of ending circle
        Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    end

    function [Cd, Pd, Cp, Pp] = decomposePath2(T) % - ... -
        % find relevant points based on T vector
        Tnorm = norm(T);
        Tunit = T/Tnorm;
        
        wp = -cross(ap, cross(Tunit, ap));
        yp = -cross(Tunit, cross(Tunit, ap));
        Cp = op + r * wp / norm(wp); % center of starting circle
        Pp = Cp - r * yp / norm(yp); % point of leaving starting circle
        
        wd = cross(ad, cross(Tunit, ad));
        yd = cross(Tunit, cross(Tunit, ad));
        Cd = od + r * wd / norm(wd); % center of ending circle
        Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    end

    function [Cd, Pd, Cp, Pp] = decomposePath3(T) % + ... -
        % find relevant points based on T vector
        Tnorm = norm(T);
        Tunit = T/Tnorm;
        
        wp = cross(ap, cross(Tunit, ap));
        yp = cross(Tunit, cross(Tunit, ap));
        Cp = op + r * wp / norm(wp); % center of starting circle
        Pp = Cp - r * yp / norm(yp); % point of leaving starting circle
        
        wd = cross(ad, cross(Tunit, ad));
        yd = cross(Tunit, cross(Tunit, ad));
        Cd = od + r * wd / norm(wd); % center of ending circle
        Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    end

    function [Cd, Pd, Cp, Pp] = decomposePath4(T) % - ... +
        % find relevant points based on T vector
        Tnorm = norm(T);
        Tunit = T/Tnorm;
        
        wp = -cross(ap, cross(Tunit, ap));
        yp = -cross(Tunit, cross(Tunit, ap));
        Cp = op + r * wp / norm(wp); % center of starting circle
        Pp = Cp - r * yp / norm(yp); % point of leaving starting circle
        
        wd = -cross(ad, cross(Tunit, ad));
        yd = -cross(Tunit, cross(Tunit, ad));
        Cd = od + r * wd / norm(wd); % center of ending circle
        Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    end

    function [theta1, gamma1, theta2, gamma2] = pathSpecs(Cd, Pd, Cp, Pp)
        W1 = Cp - op;
        w1 = W1/norm(W1);
        W2 = Cd - od;
        w2 = W2/norm(W2);
        T = Pd - Pp;
        Tunit = T / norm(T);
        theta1 = SignedAngle(ap, Tunit, cross(ap, w1));
        gamma1 = SignedAngle(bp, w1, ap);
        theta2 = SignedAngle(Tunit, ad, cross(ad, w2));
        gamma2 = SignedAngle(bd, w2, ad);
    end

% 
% -------------------------------------------------------------------

%     function err = dubinspath(T)
%         % alternate formulation based just on vector math
%         [Cd, Pd, Cp, Pp] = decomposePath(T);
%         err = T - (Pp-Pd);
%     end

    function err = dubinspath1(T) % + ... +
        % equations to solve
        Tnorm = norm(T);
        Tunit = T/Tnorm;

%         theta1 = atan2(norm(cross(ap,Tunit)),dot(ap,Tunit));
%         theta2 = atan2(norm(cross(Tunit,ad)),dot(Tunit,ad));
        theta1 = SignedAngle(ap, Tunit, cross(ap, Tunit));
        theta2 = SignedAngle(Tunit, ad, cross(Tunit, ad));

        err = T + ...
            r*(tan(theta1/2) + tan(theta2/2))*Tunit + ...
            r*(tan(theta1/2)*ap + tan(theta2/2)*ad) + ...
            - (od-op);
    end

    function err = dubinspath2(T) % - ... -
        % equations to solve
        Tnorm = norm(T);
        Tunit = T/Tnorm;

%         theta1 = atan2(norm(cross(ap,Tunit)),dot(ap,Tunit));
%         theta2 = atan2(norm(cross(Tunit,ad)),dot(Tunit,ad));
        theta1 = SignedAngle(ap, Tunit, cross(ap, Tunit));
        theta2 = SignedAngle(Tunit, ad, cross(Tunit, ad));

        err = T - ...
            r*(tan(theta1/2) + tan(theta2/2))*Tunit - ...
            r*(tan(theta1/2)*ap + tan(theta2/2)*ad) + ...
            - (od-op);
    end

    function err = dubinspath3(T) % + ... -
        % equations to solve
        Tnorm = norm(T);
        Tunit = T/Tnorm;

%         theta1 = atan2(norm(cross(ap,Tunit)),dot(ap,Tunit));
%         theta2 = atan2(norm(cross(Tunit,ad)),dot(Tunit,ad));
        theta1 = SignedAngle(ap, Tunit, cross(ap, Tunit));
        theta2 = SignedAngle(Tunit, ad, cross(Tunit, ad));

        err = T + ...
            r*(tan(theta1/2) - tan(theta2/2))*Tunit + ...
            r*(tan(theta1/2)*ap - tan(theta2/2)*ad) + ...
            - (od-op);
    end

    function err = dubinspath4(T) % - ... +
        % equations to solve
        Tnorm = norm(T);
        Tunit = T/Tnorm;

%         theta1 = atan2(norm(cross(ap,Tunit)),dot(ap,Tunit));
%         theta2 = atan2(norm(cross(Tunit,ad)),dot(Tunit,ad));
        theta1 = SignedAngle(ap, Tunit, cross(ap, Tunit));
        theta2 = SignedAngle(Tunit, ad, cross(Tunit, ad));

        err = T + ...
            r*(-tan(theta1/2) + tan(theta2/2))*Tunit + ...
            r*(-tan(theta1/2)*ap + tan(theta2/2)*ad) + ...
            - (od-op);
    end

    % Plotting without overriding
    function plotCirclePath(Op, Od, Cp, Pp, Cd, Pd, ap, ad, col1, col2)
        
        hold on
        plot3(Op(1),Op(2),Op(3), ['*' col2])
        plot3(Op(1)+[0,r*ap(1)],Op(2)+[0,r*ap(2)],Op(3)+[0,r*ap(3)],['-' col2])
        plot3(Cp(1),Cp(2),Cp(3), ['*' col2])
        plot3(Pp(1),Pp(2),Pp(3),'*g')
        
        n = cross(Op-Cp, Pp-Cp);
        v = null(n);
        theta = linspace(0,2*pi,50);
        circle_pts = bsxfun(@plus, Cp(:), r * (v(:,1)*cos(theta) + v(:,2)*sin(theta)));
        plot3(circle_pts(1,:),circle_pts(2,:),circle_pts(3,:), ['--' col2]);
        
        plot3(Od(1),Od(2),Od(3), ['*' col1])
        plot3(Od(1)+[0,r*ad(1)],Od(2)+[0,r*ad(2)],Od(3)+[0,r*ad(3)],['-' col1])
        plot3(Cd(1),Cd(2),Cd(3), ['*' col1])
        plot3(Pd(1),Pd(2),Pd(3),'*g')
        
        n = cross(Od-Cd, Pd-Cd);
        v = null(n);
        theta = linspace(0,2*pi,50);
        circle_pts = bsxfun(@plus, Cd(:), r * (v(:,1)*cos(theta) + v(:,2)*sin(theta)));
        plot3(circle_pts(1,:),circle_pts(2,:),circle_pts(3,:), ['--' col1]);
        
        
        plot3([Pd(1) Pp(1)], [Pd(2) Pp(2)], [Pd(3) Pp(3)],'-g')
    end

end