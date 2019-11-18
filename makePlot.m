points = csvread('SplinePatchToFlatPattern/points.csv');
splines = csvread('SplinePatchToFlatPattern/out.csv');

hold off
plot3(points(:,1), points(:,2), points(:,3), 'r+')
hold on;
plot3(points(:,4), points(:,5), points(:,6), 'b+')
plot3(splines(:,1), splines(:,2), splines(:,3), 'm-')
plot3(splines(:,4), splines(:,5), splines(:,6), 'c-')
axis equal
