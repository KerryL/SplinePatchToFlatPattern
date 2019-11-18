function test

points = csvread('SplinePatchToFlatPattern/points.csv');
points = points(:,1:3);
points(:,2) = abs(points(:,2));

a = round(size(points,1) / 3);
iPoints = [points(1,:); points(a+1,:); points(2 * a+1,:); points(end,:)];
%iPoints

% A general Bezier curve has four control points.  The curve intersects the
% end points, and the two inside points control the slope and curvature.
% Here, we force control points for adjacenet segments to be colinear and
% equidistant from the intersection point.

% Initial guess
cVectors = 4 * [0 1 0; points(a,:) - points(a + 2,:); points(2*a,:) - points(2*a + 2,:); 0 1 0];
%cVectors

goalFunc = @(x) doOptimization(iPoints, points, x);
x = fminsearch(goalFunc, [cVectors(1,2) cVectors(2,:) cVectors(3,:) cVectors(4,2)]);
cVectors = buildCVectors(x) %#ok<NOPRT>

drawSpline(iPoints, cVectors);
hold on
plot3(points(:,1), points(:,2), points(:,3), 'r+')
axis equal

return

function v = buildCVectors(x)
v = [0 abs(x(1)) 0; x(2:4); x(5:7); 0 abs(x(8)) 0];
return

function score = doOptimization(iPoints, goalPoints, guess)
cVectors = buildCVectors(guess);
score = computeError(iPoints, cVectors, goalPoints);
return

function p = computeSpline(iPoints, cVectors, resolution)

segments = size(iPoints,1) - 1;
p = zeros(segments * resolution, 3);
j = 1;
for i = 1:1:segments
    % TODO:  Do something smarter so we move constant increment along
    % spline, not constant increment in t
    for t = 0:1/resolution:1
        P0 = iPoints(i,:);
        if i == 1
            P1 = iPoints(i,:) + cVectors(i,:);
        else
            P1 = iPoints(i,:) - cVectors(i,:);
        end
        P2 = iPoints(i+1,:) + cVectors(i+1,:);
        P3 = iPoints(i+1,:);
        p(j,:) = (1 - t)^3 * P0 + 3 * (1 - t)^2 * t * P1 + 3 * (1 - t) * t^2 * P2 + t^3 * P3;
        j = j + 1;
    end
end

return

function drawSpline(iPoints, cVectors)

curve = computeSpline(iPoints, cVectors, 100);
hold off
plot3(curve(:,1), curve(:,2), curve(:,3), 'b-')

return

function e = computeError(iPoints, cVectors, goalPoints)

curve = computeSpline(iPoints, cVectors, 1000);
e = 0;
for i = 1:1:size(goalPoints, 1)
    xDist = (curve(:,1) - goalPoints(i,1)).^2;
    yDist = (curve(:,2) - goalPoints(i,2)).^2;
    zDist = (curve(:,3) - goalPoints(i,3)).^2;
    e = e + min(xDist + yDist + zDist);
end

return






