function [TdirMin, TmagMin, theta1min, theta2min] = dubinsCSC(r, Od, Op)
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
% Lasted edited 5/19/2023
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


% solve for optimal path specification [theta1; gamma1; theta2; gamma2]
% in each of the 4 sign combinations for gamma1, gamma2

THRESHOLD = 0.00001;

lbPP = [0, 0, 0, 0];
ubPP = [2*pi, pi, 2*pi, pi];
x0PP = [0, pi/2, 0, pi/2];
[pathSpecPP, fvalPP] = fmincon(@pathErrorNormCSC, x0PP, ...
                                [], [], [], [], lbPP, ubPP)
if abs(fvalPP) < THRESHOLD
    pathLenPP = pathLength(pathSpecPP);
else
    pathLenPP = Inf;
end

pathErrorVectorCSC(pathSpecPP)

lbPN = [0, 0, 0, -pi];
ubPN = [2*pi, pi, 2*pi, 0];
x0PN = [0, pi/2, 0, -pi/2];
[pathSpecPN, fvalPN] = fmincon(@pathErrorNormCSC, x0PN, ...
                                [], [], [], [], lbPN, ubPN)
if abs(fvalPN) < THRESHOLD
    pathLenPN = pathLength(pathSpecPN);
else
    pathLenPN = Inf;
end

pathErrorVectorCSC(pathSpecPN)


lbNP = [0, -pi, 0, 0];
ubNP = [2*pi, 0, 2*pi, pi];
x0NP = [0, pi/2, 0, pi/2];
[pathSpecNP, fvalNP] = fmincon(@pathErrorNormCSC, x0NP, ...
                                [], [], [], [], lbNP, ubNP)
if abs(fvalNP) < THRESHOLD
    pathLenNP = pathLength(pathSpecNP);
else
    pathLenNP = Inf;
end

pathErrorVectorCSC(pathSpecNP)


lbNN = [0, -pi, 0, -pi];
ubNN = [2*pi, 0, 2*pi, 0];
x0NN = [0, pi/2, 0, pi/2];
[pathSpecNN, fvalNN] = fmincon(@pathErrorNormCSC, x0NN, ...
                                [], [], [], [], lbNN, ubNN)
if abs(fvalNN) < 0.00001
    pathLenNN = pathLength(pathSpecNN);
else
    pathLenNN = Inf;
end

pathErrorVectorCSC(pathSpecNN)

pathSpecSolutions = [pathSpecPP; 
                     pathSpecPN; 
                     pathSpecNP;
                     pathSpecNN]

pathLengths = [pathLenPP;
               pathLenPN;
               pathLenNP;
               pathLenNN]

[minLen, minIdx] = min(pathLengths)


% solve for optimal t = [tx, ty, tz]
% T0 = od-op+randn(size(od));
% 
% % Initial value for Direction and Magnitude of T
% TdirectionMagnitude0 = ones(1, 4);
% TdirectionMagnitude0(1, 1:3) = ap;
% 
% % Initialize vectors for error and solutions
% TdirectionMagnitudeSolutions = zeros(4, 4);
% 
% 
% % provides individual solutions 
% TdirectionMagnitudeSolutions(1, :) = fsolve(@dubinspath1, TdirectionMagnitude0); % + ... +
% TdirectionMagnitudeSolutions(2, :) = fsolve(@dubinspath2, TdirectionMagnitude0); % - ... -
% TdirectionMagnitudeSolutions(3, :) = fsolve(@dubinspath3, TdirectionMagnitude0); % + ... -
% TdirectionMagnitudeSolutions(4, :) = fsolve(@dubinspath4, TdirectionMagnitude0); % - ... +
% 
% % double check the error is low
% SolutionError = zeros(4, 4);
% SolutionError(1, :) = dubinspath1(TdirectionMagnitudeSolutions(1, :));
% SolutionError(2, :) = dubinspath2(TdirectionMagnitudeSolutions(2, :));
% SolutionError(3, :) = dubinspath3(TdirectionMagnitudeSolutions(3, :));
% SolutionError(4, :) = dubinspath4(TdirectionMagnitudeSolutions(4, :));
% ValidPath = vecnorm(SolutionError') < 0.0001; %Threshold needs tuning?
% 
% % Display solutions and error
% disp(TdirectionMagnitudeSolutions)
% disp(SolutionError)
% 
% thetamatrix = zeros(4, 2);
% ds = zeros(4, 1);
% T_hat = zeros(4, 3);
% 
% minds = Inf;
% indexds = -1;
% 
% for i = 1:4
% 
%     Tdir = TdirectionMagnitudeSolutions(i, 1:3);
%     Tmag = TdirectionMagnitudeSolutions(i, 4);
%     T_hat(i, :) = Tdir / norm(Tdir);
%     T_found = Tmag * T_hat(i, :);
% 
%     thetamatrix(i, 1) = SignedAngle(ap, T_hat(i, :), cross(ap, T_hat(i, :)));
% 
%     % find the point om after a C + S
%     wp = cross(ap, T_hat(i, :));
%     if ap == T_hat(i, :)
%         % No turning in first segment
%         om = op' + T_found';
%     else
%         % rotate ap about the circle, then move forward by T_found
%         wp = wp/norm(wp);
%         o2c1 = -cross(ap, wp)';
%         o2c1 = o2c1/norm(o2c1);
%         om = op' + r*(eye(3) - RotationalMatrix(wp, thetamatrix(i, 1)))* o2c1 + T_found';
%     end
%     % find the vector from om to od, create a normal vector as the axis of
%     % rotation
%     if ValidPath(i)
%         om2od = od' - om;
%         ntemp = cross(T_hat(i, :), om2od);
% 
%         thetamatrix(i, 2) = SignedAngle(T_hat(i, :), ad, ntemp);
% 
%         ds(i) = r*abs(thetamatrix(i, 1)) + Tmag + r*abs(thetamatrix(i, 2));
%         if ds(i) < minds
%             minds = ds(i);
%             indexds = i;
%         end
%     else
%         ds(i) = Inf;
%     end
% 
% end
% 
% disp(minds);
% 
% % Output minimum Tunittor as well as minimum theta1 and theta2 values
% TdirMin = TdirectionMagnitudeSolutions(indexds, 1:3);
% TmagMin = TdirectionMagnitudeSolutions(indexds, 4);
% theta1min = thetamatrix(indexds, 1);
% theta2min = thetamatrix(indexds, 2);
% 
% disp(thetamatrix)

% Plotting
% --------------------------------------------------------------------

% decompose path in relevant points for plotting
%[Cd1, Pd1, Cp1, Pp1] = decomposePath1(TdirectionMagnitudeSolutions(1,:));
[CdPP, PdPP, tdPP, CpPP, PpPP, tpPP] = constructPathCSC(pathSpecPP)
[CdPN, PdPN, tdPN, CpPN, PpPN, tpPN] = constructPathCSC(pathSpecPN)
[CdNP, PdNP, tdNP, CpNP, PpNP, tpNP] = constructPathCSC(pathSpecNP)
[CdNN, PdNN, tdNN, CpNN, PpNN, tpNN] = constructPathCSC(pathSpecNN)

% plotting
hold on
subplot(2, 2, 1)
%plotCirclePath(op, od, Cp1, Pp1, Cd1, Pd1, ap, ad, 'b', 'r')
plotCirclePath(op, od, CpPP, PpPP, CdPP, PdPP, ap, ad, 'c', 'm');

axis equal

% decompose path in relevant points for plotting
%[Cd2, Pd2, Cp2, Pp2] = decomposePath2(TdirectionMagnitudeSolutions(2,:));


% plotting
hold on
subplot(2, 2, 2)
%plotCirclePath(op, od, Cp2, Pp2, Cd2, Pd2, ap, ad, 'b', 'r')
plotCirclePath(op, od, CpNN, PpNN, CdNN, PdNN, ap, ad, 'c', 'm');


axis equal

% decompose path in relevant points for plotting
%[Cd3, Pd3, Cp3, Pp3] = decomposePath3(TdirectionMagnitudeSolutions(3,:));


% plotting
hold on
subplot(2, 2, 3)
%plotCirclePath(op, od, Cp3, Pp3, Cd3, Pd3, ap, ad, 'b', 'r')
plotCirclePath(op, od, CpPN, PpPN, CdPN, PdPN, ap, ad, 'c', 'm');

axis equal

% decompose path in relevant points for plotting
%[Cd4, Pd4, Cp4, Pp4] = decomposePath4(TdirectionMagnitudeSolutions(4,:));

% plotting
hold on
subplot(2, 2, 4)
%plotCirclePath(op, od, Cp4, Pp4, Cd4, Pd4, ap, ad, 'b', 'r')
plotCirclePath(op, od, CpNP, PpNP, CdNP, PdNP, ap, ad, 'c', 'm');


axis equal

AllPathSpecs = zeros(4,4);
[th11,g11,th21,g21] = pathSpecs(CdPP, PdPP, CpPP, PpPP);
AllPathSpecs(1,:) = [th11,g11,th21,g21];
[th12,g12,th22,g22] = pathSpecs(CdNN, PdNN, CpNN, PpNN);
AllPathSpecs(2,:) = [th12,g12,th22,g22];
[th13,g13,th23,g23] = pathSpecs(CdPN, PdPN, CpPN, PpPN);
AllPathSpecs(3,:) = [th13,g13,th23,g23];
[th14,g14,th24,g24] = pathSpecs(CdNP, PdNP, CpNP, PpNP);
AllPathSpecs(4,:) = [th14,g14,th24,g24];
AllPathSpecs

% Functions
% ---------------------------------------------------------------------

    % function [Cd, Pd, Cp, Pp] = decomposePath1(TdirectionMagnitude) % + ... +
    %     % find relevant points based on T vector
    %     Tdirection = TdirectionMagnitude(1:3);
    %     Tunit = Tdirection / norm(Tdirection);
    % 
    %     wp = cross(ap, cross(Tunit, ap));
    %     yp = cross(Tunit, cross(Tunit, ap));
    %     Cp = op + r * wp / norm(wp); % center of starting circle
    %     Pp = Cp - r * yp / norm(yp); % point of leaving starting circle
    % 
    %     wd = -cross(ad, cross(Tunit, ad));
    %     yd = -cross(Tunit, cross(Tunit, ad));
    %     Cd = od + r * wd / norm(wd); % center of ending circle
    %     Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    % end
    % 
    % function [Cd, Pd, Cp, Pp] = decomposePath2(TdirectionMagnitude) % - ... -
    %     % find relevant points based on T vector
    %     Tdirection = TdirectionMagnitude(1:3);
    %     Tunit = Tdirection / norm(Tdirection);
    % 
    %     wp = -cross(ap, cross(Tunit, ap));
    %     yp = -cross(Tunit, cross(Tunit, ap));
    %     Cp = op + r * wp / norm(wp); % center of starting circle
    %     Pp = Cp - r * yp / norm(yp); % point of leaving starting circle
    % 
    %     wd = cross(ad, cross(Tunit, ad));
    %     yd = cross(Tunit, cross(Tunit, ad));
    %     Cd = od + r * wd / norm(wd); % center of ending circle
    %     Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    % end
    % 
    % function [Cd, Pd, Cp, Pp] = decomposePath3(TdirectionMagnitude) % + ... -
    %     % find relevant points based on T vector
    %     Tdirection = TdirectionMagnitude(1:3);
    %     Tunit = Tdirection / norm(Tdirection);
    % 
    %     wp = cross(ap, cross(Tunit, ap));
    %     yp = cross(Tunit, cross(Tunit, ap));
    %     Cp = op + r * wp / norm(wp); % center of starting circle
    %     Pp = Cp - r * yp / norm(yp); % point of leaving starting circle
    % 
    %     wd = cross(ad, cross(Tunit, ad));
    %     yd = cross(Tunit, cross(Tunit, ad));
    %     Cd = od + r * wd / norm(wd); % center of ending circle
    %     Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    % end
    % 
    % function [Cd, Pd, Cp, Pp] = decomposePath4(TdirectionMagnitude) % - ... +
    %     % find relevant points based on T vector
    %     Tdirection = TdirectionMagnitude(1:3);
    %     Tunit = Tdirection / norm(Tdirection);
    % 
    %     wp = -cross(ap, cross(Tunit, ap));
    %     yp = -cross(Tunit, cross(Tunit, ap));
    %     Cp = op + r * wp / norm(wp); % center of starting circle
    %     Pp = Cp - r * yp / norm(yp); % point of leaving starting circle
    % 
    %     wd = -cross(ad, cross(Tunit, ad));
    %     yd = -cross(Tunit, cross(Tunit, ad));
    %     Cd = od + r * wd / norm(wd); % center of ending circle
    %     Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    % end

% 
% -------------------------------------------------------------------

%     function err = dubinspath(T)
%         % alternate formulation based just on vector math
%         [Cd, Pd, Cp, Pp] = decomposePath(T);
%         err = T - (Pp-Pd);
%     end

    % Given the path specifications [theta1; gamma1; theta2; gamma2] with:
    %   Path arc angles theta1, theta2 each in [0,2pi]
    %   Circle orientation angles gamma1, gamma2 each in [-pi,pi]
    %       measured about ap, ad from bp, bd respectively
    %       i.e., if you rotate bp about ap by gamma1, the result will
    %       point from op towards Cp (the circle 1 center)
    %       and similarly for gamma2
    % Output:
    %   Cp, Cd      Circle centers
    %   Pp, Pd      Points where the S component meets the circles
    %   tp, td      Proximal and distal unit tangents 
    %               Proximal is the rotation of op by theta1 about circle 1
    %               Distal is rotation of od by -theta1 about circle 2 
    %               (negative because theta2 is measured from Pd to od)
    %               Note that in a valid path the tangents align tp == td

    function [Cd, Pd, td, Cp, Pp, tp] = constructPathCSC(pathSpecification)
        theta1 = pathSpecification(1);
        gamma1 = pathSpecification(2);
        theta2 = pathSpecification(3);
        gamma2 = pathSpecification(4);
        w1 = (RotationalMatrix(ap, gamma1) * bp')'; % direction from op to Cp
        w1 = w1 / norm(w1);
        Cp = op + r * w1;
        w2 = (RotationalMatrix(ad, gamma2) * bd')'; % direction from od to Cd
        w2 = w2 / norm(w2);
        Cd = od + r * w2;

        % Note that w1 faces inwards while ap faces forwards.
        % By the right-hand rules for cross product and curve orientation, 
        % cross(ap, w1) about Cp is normal to circle 1 such that
        % ap aligns with direction of increasing theta.
        % Similarly, measuring angles about cross(ad, w2) aligns with ad.

        R1 = RotationalMatrix(cross(ap, w1), theta1); % Rotation along circle 1
        R2inv = RotationalMatrix(cross(ad, w2), -theta2); % Backwards along circle 2
        
        tp = (R1 * ap')';        
        y1 = (R1 * w1')'; %rotate from w1
        Pp = Cp - r*y1;

        td = (R2inv * ad')';
        y2 = (R2inv * w2')'; %rotate from w2
        Pd = Cd - r*y2;
    end

    % Given path specification vector [theta1; gamma1; theta2; gamma2],
    % return the length of the corresponding path
    function len = pathLength(pathSpecification)
        theta1 = pathSpecification(1);
        theta2 = pathSpecification(3);
        [Cd, Pd, td, Cp, Pp, tp] = constructPathCSC(pathSpecification);
        len = r*theta1 + norm(Pd-Pp) + r*theta2;
    end
    

    % pathSpecification is a vector [theta1; gamma1; theta2; gamma2] with:
    %   Path arc angles theta1, theta2 each in [0,2pi)
    %   Circle orientation angles gamma1, gamma2
    %       measured about ap, ad from bp, bd respectively
    %       i.e., if you rotate bp about ap by gamma1, the result will
    %       point from op towards Cp (the circle 1 center)
    %       and similarly for gamma2
    %       This objective will be used in fmincon with constraint ranges 
    %       on gamma1, gamma2 varying to cover the 4 path type cases
    % The error function measures:
    %   Alignment between tangents to the circles   (3 equations)
    %   Point-to-line displacement                  (3 equations)
    %       from Pd to the line through Pp in direction tp

    function err = pathErrorVectorCSC(pathSpecification)
        [Cd, Pd, td, Cp, Pp, tp] = constructPathCSC(pathSpecification);
        err = zeros(6,1);
        err(1:3) = td - tp; % Tangent Alignment Error
        err(4:6) = Pd - nearestPointOnRay(Pd, Pp, tp); % Point Alignment Error        
    end

    % fmincon expects f to output a scalar that it will try to minimize
    % since the error vector should be [0,0,0,0,0,0] in a valid solution,
    % we can return the error vector's norm as the objective to minimize
    % and consider a solution value iff it gets the objective to 0
    function obj = pathErrorNormCSC(pathSpecification)
        obj = norm(pathErrorVectorCSC(pathSpecification));        
    end

    % Implement the constraints 
    %       gamma1 >= min1      i.e.    -gamma1 <= -min1
    %       gamma1 <= max1
    %       gamma2 >= min2      i.e.    -gamma2 <= -min2
    %       gamma2 <= max2
    % as matrix A and vector b for constraint equation Ax <= b
    %       where x = [theta1; gamma1; theta2; gamma2]
    % function [A, b] = gammaRanges(min1, max1, min2, max2)
    %     A = [0,  -1,   0,   0;
    %          0,   1,   0,   0;
    %          0,   0,   0,  -1;
    %          0,   0,   0,   1];
    %     b = [-min1; 
    %           max1; 
    %          -min2; 
    %           max2];
    % end
   
    
    % function err = dubinspath1(TdirectionMagnitude) % + ... +
    %     % equations to solve
    %     Tdir = TdirectionMagnitude(1:3);
    %     Tmag = TdirectionMagnitude(4);
    %     Tunit = Tdir / norm(Tdir);
    % 
    %     theta1 = SignedAngle(ap, Tunit, cross(ap, Tunit));
    %     theta2 = SignedAngle(Tunit, ad, cross(Tunit, ad));
    % 
    %     T_err = Tmag * Tunit + ...
    %         r*(tan(theta1/2) + tan(theta2/2))*Tunit + ...
    %         r*(tan(theta1/2)*ap + tan(theta2/2)*ad) + ...
    %         - (od-op);
    % 
    %     err = zeros(1,4);
    %     err(1:3) = T_err;
    %     err(4) = norm(Tdir) - 1; % we want Tdir == Tunit
    % end
    % 
    % function err = dubinspath2(TdirectionMagnitude) % + ... +
    %     % equations to solve
    %     Tdir = TdirectionMagnitude(1:3);
    %     Tmag = TdirectionMagnitude(4);
    %     Tunit = Tdir / norm(Tdir);
    % 
    %     theta1 = SignedAngle(ap, Tunit, cross(ap, Tunit));
    %     theta2 = SignedAngle(Tunit, ad, cross(Tunit, ad));
    % 
    %     T_err = Tmag * Tunit - ...
    %         r*(tan(theta1/2) + tan(theta2/2))*Tunit - ...
    %         r*(tan(theta1/2)*ap + tan(theta2/2)*ad) + ...
    %         - (od-op);
    % 
    %     err = zeros(1,4);
    %     err(1:3) = T_err;
    %     err(4) = norm(Tdir) - 1; % we want Tdir == Tunit
    % end
    % 
    % function err = dubinspath3(TdirectionMagnitude) % + ... +
    %     % equations to solve
    %     Tdir = TdirectionMagnitude(1:3);
    %     Tmag = TdirectionMagnitude(4);
    %     Tunit = Tdir / norm(Tdir);
    % 
    %     theta1 = SignedAngle(ap, Tunit, cross(ap, Tunit));
    %     theta2 = SignedAngle(Tunit, ad, cross(Tunit, ad));
    % 
    %     T_err = Tmag * Tunit + ...
    %         r*(tan(theta1/2) - tan(theta2/2))*Tunit + ...
    %         r*(tan(theta1/2)*ap - tan(theta2/2)*ad) + ...
    %         - (od-op);
    % 
    %     err = zeros(1,4);
    %     err(1:3) = T_err;
    %     err(4) = norm(Tdir) - 1; % we want Tdir == Tunit
    % end
    % 
    % function err = dubinspath4(TdirectionMagnitude) % + ... +
    %     % equations to solve
    %     Tdir = TdirectionMagnitude(1:3);
    %     Tmag = TdirectionMagnitude(4);
    %     Tunit = Tdir / norm(Tdir);
    % 
    %     theta1 = SignedAngle(ap, Tunit, cross(ap, Tunit));
    %     theta2 = SignedAngle(Tunit, ad, cross(Tunit, ad));
    % 
    %     T_err = Tmag * Tunit + ...
    %         r*(-tan(theta1/2) + tan(theta2/2))*Tunit + ...
    %         r*(-tan(theta1/2)*ap + tan(theta2/2)*ad) + ...
    %         - (od-op);
    % 
    %     err = zeros(1,4);
    %     err(1:3) = T_err;
    %     err(4) = norm(Tdir) - 1; % we want Tdir == Tunit
    % end

   % For debugging
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