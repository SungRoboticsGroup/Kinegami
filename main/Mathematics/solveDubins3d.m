function [TdirMin, TmagMin, theta1min, theta2min] = solveDubins3d(r, Od, Op)
% SOLVEDUBINS3D - Find the CSC Dubins path from frame Op to Od.
% Outputs row vector and two radian values. Ensure that Od and Op are
% normalized prior to entry. Only guaranteed to work if op and od are
% at least 4r apart. Based on S. Hota and D. Ghose, 
% "Optimal Geometrical Path in 3D with Curvature Constraint", IROS 2010.

% Inputs:
%   r           - desired radius of folded origami linkage.
%   Od          - distal, or terminal, frame of dubins pathing.
%   Op          - proximal, or initial, frame of dubins pathing.

% Outputs:
%   TdirMin     - unit direction vector of the S part of the shortest CSC 
%                 path for Dubins linkage. If the S part is empty, this
%                 represents the common tangent connecting the circles.
%   TmagMin     - magnitude of the S part of the shortest CSC path.
%   theta1min   - theta value (in radians) for use in RotationalMatrix.m,
%                 for utility in determining rotational matrix around wp.
%   theta2min   - theta value (in radians) for use in RotationalMatrix.m,
%                 for utility in determining rotational matrix around wd.

% Authors: 
% Cynthia Sung <crsung@seas.upenn.edu>
% Lucien Peach <peach@seas.upenn.edu>
% Wei-Hsi Chen <weicc@seas.upenn.edu>
% Daniel Feshbach <feshbach@seas.upenn.edu>
% Lasted edited 6/5/2023
%
% Copyright (C) 2023 The Trustees of the University of Pennsylvania. 
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
T0 = od - op + r*ap; % randn(size(od));
if norm(T0)==0
    T0 = r*ap;
end
% 
% Initial value for Direction and Magnitude of T
TdirectionMagnitude0 = ones(1, 4);
TdirectionMagnitude0(1, 1:3) = T0 / norm(T0);
TdirectionMagnitude0(1, 4) = norm(T0);

% Initialize vectors for error and solutions
TdirectionMagnitudeSolutions = zeros(4, 4);


% provides individual solutions 
TdirectionMagnitudeSolutions(1, :) = fsolve(@dubinspath1, TdirectionMagnitude0); % + ... +
TdirectionMagnitudeSolutions(2, :) = fsolve(@dubinspath2, TdirectionMagnitude0); % - ... -
TdirectionMagnitudeSolutions(3, :) = fsolve(@dubinspath3, TdirectionMagnitude0); % + ... -
TdirectionMagnitudeSolutions(4, :) = fsolve(@dubinspath4, TdirectionMagnitude0); % - ... +

% Absolutize Tmag
TdirectionMagnitudeSolutions(:,4) = abs(TdirectionMagnitudeSolutions(:,4));

% double check the error is low
SolutionError = zeros(4, 4);
SolutionError(1, :) = dubinspath1(TdirectionMagnitudeSolutions(1, :));
SolutionError(2, :) = dubinspath2(TdirectionMagnitudeSolutions(2, :));
SolutionError(3, :) = dubinspath3(TdirectionMagnitudeSolutions(3, :));
SolutionError(4, :) = dubinspath4(TdirectionMagnitudeSolutions(4, :));
ValidPath = vecnorm(SolutionError') < 0.0001; %Threshold needs tuning?

% Display solutions and error
disp(TdirectionMagnitudeSolutions)
disp(SolutionError)

thetamatrix = zeros(4, 2);
ds = zeros(4, 1);
T_hat = zeros(4, 3);

minds = Inf;
indexds = -1;

for i = 1:4

    Tdir = TdirectionMagnitudeSolutions(i, 1:3);
    Tmag = TdirectionMagnitudeSolutions(i, 4);
    T_hat(i, :) = Tdir / norm(Tdir);
    T_found = Tmag * T_hat(i, :);

    thetamatrix(i, 1) = SignedAngle(ap, T_hat(i, :), cross(ap, T_hat(i, :)));

    % find the point om after a C + S
    wp = cross(ap, T_hat(i, :));
    if ap == T_hat(i, :)
        % No turning in first segment
        om = op' + T_found';
    else
        % rotate ap about the circle, then move forward by T_found
        wp = wp/norm(wp);
        o2c1 = -cross(ap, wp)';
        o2c1 = o2c1/norm(o2c1);
        om = op' + r*(eye(3) - RotationalMatrix(wp, thetamatrix(i, 1)))* o2c1 + T_found';
    end
    % find the vector from om to od, create a normal vector as the axis of
    % rotation
    if ValidPath(i)
        om2od = od' - om;
        ntemp = cross(T_hat(i, :), om2od);

        thetamatrix(i, 2) = SignedAngle(T_hat(i, :), ad, ntemp);

        ds(i) = r*abs(thetamatrix(i, 1)) + Tmag + r*abs(thetamatrix(i, 2));
        if ds(i) < minds
            minds = ds(i);
            indexds = i;
        end
    else
        ds(i) = Inf;
    end

end

disp(minds);

% Output minimum Tunittor as well as minimum theta1 and theta2 values
TdirMin = TdirectionMagnitudeSolutions(indexds, 1:3);
TmagMin = TdirectionMagnitudeSolutions(indexds, 4);
theta1min = thetamatrix(indexds, 1);
theta2min = thetamatrix(indexds, 2);

disp(thetamatrix)

% Plotting
% --------------------------------------------------------------------

% decompose path in relevant points for plotting
[Cd1, Pd1, Cp1, Pp1] = decomposePath1(TdirectionMagnitudeSolutions(1,:));


% plotting
hold on
subplot(2, 2, 1)
plotCirclePath(op, od, Cp1, Pp1, Cd1, Pd1, ap, ad, 'b', 'r')

axis equal

% decompose path in relevant points for plotting
[Cd2, Pd2, Cp2, Pp2] = decomposePath2(TdirectionMagnitudeSolutions(2,:));


% plotting
hold on
subplot(2, 2, 2)
plotCirclePath(op, od, Cp2, Pp2, Cd2, Pd2, ap, ad, 'b', 'r')

axis equal

% decompose path in relevant points for plotting
[Cd3, Pd3, Cp3, Pp3] = decomposePath3(TdirectionMagnitudeSolutions(3,:));


% plotting
hold on
subplot(2, 2, 3)
plotCirclePath(op, od, Cp3, Pp3, Cd3, Pd3, ap, ad, 'b', 'r')

axis equal

% decompose path in relevant points for plotting
[Cd4, Pd4, Cp4, Pp4] = decomposePath4(TdirectionMagnitudeSolutions(4,:));

% plotting
hold on
subplot(2, 2, 4)
plotCirclePath(op, od, Cp4, Pp4, Cd4, Pd4, ap, ad, 'b', 'r')

axis equal

% Functions
% ---------------------------------------------------------------------
    
   
    function [Cd, Pd, Cp, Pp] = decomposePath1(TdirectionMagnitude) % + ... +
        % find relevant points based on T vector
        Tdirection = TdirectionMagnitude(1:3);
        Tunit = Tdirection / norm(Tdirection);

        wp = cross(ap, unitNormalToBoth(Tunit, ap));
        yp = cross(Tunit, unitNormalToBoth(Tunit, ap));
        
        Cp = op + r * wp / norm(wp); % center of starting circle
        Pp = Cp - r * yp / norm(yp); % point of leaving starting circle
        
        wd = -cross(ad, unitNormalToBoth(Tunit, ad));
        yd = -cross(Tunit, unitNormalToBoth(Tunit, ad));
        Cd = od + r * wd / norm(wd); % center of ending circle
        Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    end

    function [Cd, Pd, Cp, Pp] = decomposePath2(TdirectionMagnitude) % - ... -
        % find relevant points based on T vector
        Tdirection = TdirectionMagnitude(1:3);
        Tunit = Tdirection / norm(Tdirection);

        wp = -cross(ap, unitNormalToBoth(Tunit, ap));
        yp = -cross(Tunit, unitNormalToBoth(Tunit, ap));
        Cp = op + r * wp / norm(wp); % center of starting circle
        Pp = Cp - r * yp / norm(yp); % point of leaving starting circle

        wd = cross(ad, unitNormalToBoth(Tunit, ad));
        yd = cross(Tunit, unitNormalToBoth(Tunit, ad));
        Cd = od + r * wd / norm(wd); % center of ending circle
        Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    end

    function [Cd, Pd, Cp, Pp] = decomposePath3(TdirectionMagnitude) % + ... -
        % find relevant points based on T vector
        Tdirection = TdirectionMagnitude(1:3);
        Tunit = Tdirection / norm(Tdirection);

        wp = cross(ap, unitNormalToBoth(Tunit, ap));
        yp = cross(Tunit, unitNormalToBoth(Tunit, ap));
        Cp = op + r * wp / norm(wp); % center of starting circle
        Pp = Cp - r * yp / norm(yp); % point of leaving starting circle

        wd = cross(ad, unitNormalToBoth(Tunit, ad));
        yd = cross(Tunit, unitNormalToBoth(Tunit, ad));
        Cd = od + r * wd / norm(wd); % center of ending circle
        Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    end

    function [Cd, Pd, Cp, Pp] = decomposePath4(TdirectionMagnitude) % - ... +
        % find relevant points based on T vector
        Tdirection = TdirectionMagnitude(1:3);
        Tunit = Tdirection / norm(Tdirection);

        wp = -cross(ap, unitNormalToBoth(Tunit, ap));
        yp = -cross(Tunit, unitNormalToBoth(Tunit, ap));
        Cp = op + r * wp / norm(wp); % center of starting circle
        Pp = Cp - r * yp / norm(yp); % point of leaving starting circle

        wd = -cross(ad, unitNormalToBoth(Tunit, ad));
        yd = -cross(Tunit, unitNormalToBoth(Tunit, ad));
        Cd = od + r * wd / norm(wd); % center of ending circle
        Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    end

% 
% -------------------------------------------------------------------

%     function err = dubinspath(T)
%         % alternate formulation based just on vector math
%         [Cd, Pd, Cp, Pp] = decomposePath(T);
%         err = T - (Pp-Pd);
%     end
    

    function err = dubinspath1(TdirectionMagnitude) % + ... +
        % equations to solve
        Tdir = TdirectionMagnitude(1:3);
        Tmag = abs(TdirectionMagnitude(4));
        Tunit = Tdir / norm(Tdir);

        theta1 = SignedAngle(ap, Tunit, cross(ap, Tunit));
        theta2 = SignedAngle(Tunit, ad, cross(Tunit, ad));

        T_err = Tmag * Tunit + ...
            r*(tan(theta1/2) + tan(theta2/2))*Tunit + ...
            r*(tan(theta1/2)*ap + tan(theta2/2)*ad) + ...
            - (od-op);

        err = zeros(1,4);
        err(1:3) = T_err;
        err(4) = norm(Tdir) - 1; % we want Tdir == Tunit
    end

    function err = dubinspath2(TdirectionMagnitude) % + ... +
        % equations to solve
        Tdir = TdirectionMagnitude(1:3);
        Tmag = abs(TdirectionMagnitude(4));
        Tunit = Tdir / norm(Tdir);

        theta1 = SignedAngle(ap, Tunit, cross(ap, Tunit));
        theta2 = SignedAngle(Tunit, ad, cross(Tunit, ad));

        T_err = Tmag * Tunit - ...
            r*(tan(theta1/2) + tan(theta2/2))*Tunit - ...
            r*(tan(theta1/2)*ap + tan(theta2/2)*ad) + ...
            - (od-op);

        err = zeros(1,4);
        err(1:3) = T_err;
        err(4) = norm(Tdir) - 1; % we want Tdir == Tunit
    end

    function err = dubinspath3(TdirectionMagnitude) % + ... +
        % equations to solve
        Tdir = TdirectionMagnitude(1:3);
        Tmag = abs(TdirectionMagnitude(4));
        Tunit = Tdir / norm(Tdir);

        theta1 = SignedAngle(ap, Tunit, cross(ap, Tunit));
        theta2 = SignedAngle(Tunit, ad, cross(Tunit, ad));

        T_err = Tmag * Tunit + ...
            r*(tan(theta1/2) - tan(theta2/2))*Tunit + ...
            r*(tan(theta1/2)*ap - tan(theta2/2)*ad) + ...
            - (od-op);

        err = zeros(1,4);
        err(1:3) = T_err;
        err(4) = norm(Tdir) - 1; % we want Tdir == Tunit
    end

    function err = dubinspath4(TdirectionMagnitude) % + ... +
        % equations to solve
        Tdir = TdirectionMagnitude(1:3);
        Tmag = abs(TdirectionMagnitude(4));
        Tunit = Tdir / norm(Tdir);

        theta1 = SignedAngle(ap, Tunit, cross(ap, Tunit));
        theta2 = SignedAngle(Tunit, ad, cross(Tunit, ad));

        T_err = Tmag * Tunit + ...
            r*(-tan(theta1/2) + tan(theta2/2))*Tunit + ...
            r*(-tan(theta1/2)*ap + tan(theta2/2)*ad) + ...
            - (od-op);

        err = zeros(1,4);
        err(1:3) = T_err;
        err(4) = norm(Tdir) - 1; % we want Tdir == Tunit
    end

   % % For debugging
   %  function [theta1, gamma1, theta2, gamma2] = pathAngleSpecs(Cd, Pd, Cp, Pp)
   %      W1 = Cp - op;
   %      w1 = W1/norm(W1);
   %      W2 = Cd - od;
   %      w2 = W2/norm(W2);
   %      T = Pd - Pp;
   %      Tunit = T / norm(T);
   %      theta1 = SignedAngle(ap, Tunit, cross(ap, w1));
   %      gamma1 = SignedAngle(bp, w1, ap);
   %      theta2 = SignedAngle(Tunit, ad, cross(ad, w2));
   %      gamma2 = SignedAngle(bd, w2, ad);
   %  end


    % Plotting without overriding
    function plotCirclePath(Op, Od, Cp, Pp, Cd, Pd, ap, ad, col1, col2)
        
        hold on
        plot3(Op(1),Op(2),Op(3), ['*' col2])
        plot3(Op(1)+[0,r*ap(1)],Op(2)+[0,r*ap(2)],Op(3)+[0,r*ap(3)],['-' col2])
        plot3(Cp(1),Cp(2),Cp(3), ['*' col2])
        plot3(Pp(1),Pp(2),Pp(3),'*g')
        
        %n = cross(Op-Cp, Pp-Cp);
        n = cross (Op-Cp, r*ap);
        v = null(n);
        theta = linspace(0,2*pi,50);
        circle_pts = bsxfun(@plus, Cp(:), r * (v(:,1)*cos(theta) + v(:,2)*sin(theta)));
        plot3(circle_pts(1,:),circle_pts(2,:),circle_pts(3,:), ['--' col2]);
        
        plot3(Od(1),Od(2),Od(3), ['*' col1])
        plot3(Od(1)+[0,r*ad(1)],Od(2)+[0,r*ad(2)],Od(3)+[0,r*ad(3)],['-' col1])
        plot3(Cd(1),Cd(2),Cd(3), ['*' col1])
        plot3(Pd(1),Pd(2),Pd(3),'*g')
        
        %n = cross(Od-Cd, Pd-Cd);
        n = cross (Od-Cd, r*ad);
        v = null(n);
        theta = linspace(0,2*pi,50);
        circle_pts = bsxfun(@plus, Cd(:), r * (v(:,1)*cos(theta) + v(:,2)*sin(theta)));
        plot3(circle_pts(1,:),circle_pts(2,:),circle_pts(3,:), ['--' col1]);
        
        
        plot3([Pd(1) Pp(1)], [Pd(2) Pp(2)], [Pd(3) Pp(3)],'-g')
    end

end