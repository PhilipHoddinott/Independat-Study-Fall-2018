clear all; close all;

% Define the input grid
[x, y] = meshgrid(linspace(-1, 1));
% Calculate the two surfaces
z1 = y.^2 + 2*x;
z2 = 2*y.^3 - x.^2;
% Visualize the two surfaces
surface(x, y, z1, 'FaceColor', [0.5 1.0 0.5], 'EdgeColor', 'none');
surface(x, y, z2, 'FaceColor', [1.0 0.5 0.0], 'EdgeColor', 'none');
view(3); camlight; axis vis3d
% Take the difference between the two surface heights and find the contour
% where that surface is zero.
zdiff = z1 - z2;
C = contours(x, y, zdiff, [0 0]);
% Extract the x- and y-locations from the contour matrix C.
xL = C(1, 2:end);
yL = C(2, 2:end);
% Interpolate on the first surface to find z-locations for the intersection
% line.
zL = interp2(x, y, z1, xL, yL);
% Visualize the line.
line(xL, yL, zL, 'Color', 'k', 'LineWidth', 3);