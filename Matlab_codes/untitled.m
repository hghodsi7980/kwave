clear; clc; close all;

% Parameters for the test case
a = 75;
b = 100;
c = 50;
numPoints = 5000;

figure;

% ==========================================================
% Method 1: Voronoi-Based Structure
% ==========================================================
subplot(2, 2, 1);
% Generate random seed points for Voronoi tessellation
numCells = randi([5, 20]); % Number of Voronoi cells
seedPoints = [a * rand(numCells, 1), b * rand(numCells, 1), c * rand(numCells, 1)];

% Compute the Voronoi tessellation
[V, C] = voronoin(seedPoints);

% Randomly select a few cells and create a convex hull for each
Points = [];
for i = 1:length(C)
    cellPoints = V(C{i}, :);
    
    % Keep only points within the bounding box (exclude infinite points)
    finitePoints = cellPoints(all(cellPoints < [a b c] & cellPoints > 0, 2), :);
    if ~isempty(finitePoints)
        Points = [Points; finitePoints];
    end
end

% Apply some random displacements to the points
if ~isempty(Points)
    displacements = [randn(size(Points, 1), 1) * a / 10, randn(size(Points, 1), 1) * b / 10, randn(size(Points, 1), 1) * c / 10];
    Points = Points + displacements;
    DT = delaunayTriangulation(Points);
    trisurf(DT.ConnectivityList, Points(:,1), Points(:,2), Points(:,3), 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
else
    disp('No points generated for Voronoi structure');
end
title('Voronoi-Based Structure');
axis equal;
view(3);

% ==========================================================
% Method 2: Perlin Noise Displacement
% ==========================================================
subplot(2, 2, 2);
% Generate the initial ellipsoid-like points
theta = 2 * pi * rand(numPoints, 1);
phi = pi * rand(numPoints, 1);
x_ellipsoid = a * sin(phi) .* cos(theta);
y_ellipsoid = b * sin(phi) .* sin(theta);
z_ellipsoid = c * cos(phi);
Points = [x_ellipsoid, y_ellipsoid, z_ellipsoid];

% Apply random displacement as a placeholder for Perlin noise
noise_scale = 5; % Scale factor for noise effect
for i = 1:size(Points, 1)
    noise_value = rand() - 0.5; % Placeholder for Perlin noise value
    Points(i, :) = Points(i, :) + noise_scale * noise_value * [1, 1, 1];
end

% Create and plot the triangulated surface
DT = delaunayTriangulation(Points);
trisurf(DT.ConnectivityList, Points(:,1), Points(:,2), Points(:,3), 'FaceColor', 'magenta', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
title('Perlin Noise Displacement');
axis equal;
view(3);

% ==========================================================
% Method 3: Fractal-Based Approach
% ==========================================================
subplot(2, 2, 3);
% Generate a grid of points
grid_size = 20;
[X, Y, Z] = meshgrid(linspace(-a, a, grid_size), linspace(-b, b, grid_size), linspace(-c, c, grid_size));
Points = [X(:), Y(:), Z(:)];

% Apply a fractal-like displacement
max_displacement = 10;
displacements = max_displacement * (rand(size(Points, 1), 3) - 0.5);
Points = Points + displacements;

% Create and plot the triangulated surface
DT = delaunayTriangulation(Points);
trisurf(DT.ConnectivityList, Points(:,1), Points(:,2), Points(:,3), 'FaceColor', 'yellow', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
title('Fractal-Based Approach');
axis equal;
view(3);

% ==========================================================
% Method 4: Random Cluster Formation
% ==========================================================
subplot(2, 2, 4);
% Number of clusters
numClusters = randi([3, 10]);
Points = [];

for i = 1:numClusters
    cluster_center = [a * rand(), b * rand(), c * rand()];
    cluster_points = mvnrnd(cluster_center, eye(3) * 50, round(numPoints / numClusters));
    Points = [Points; cluster_points];
end

% Create and plot the triangulated surface
DT = delaunayTriangulation(Points);
trisurf(DT.ConnectivityList, Points(:,1), Points(:,2), Points(:,3), 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
title('Random Cluster Formation');
axis equal;
view(3);
