% Number of total points
num_points = 1000;

% Number of dense clusters
num_clusters = 50;

% Generate random coordinates for cluster centers
cluster_centers = rand(num_clusters, 3);

% Define bandwidth for Gaussian kernel (controls the spread)
bandwidth = 0.1; % adjust as needed

% Generate points using KDE around each cluster center
points = [];
for i = 1:num_clusters
    % Generate points around each cluster center using Gaussian kernel
    cluster_points = mvnrnd(cluster_centers(i, :), bandwidth^2 * eye(3), num_points);
    points = [points; cluster_points];
end

% Plot the points
scatter3(points(:,1), points(:,2), points(:,3), '.');

% Adjust axis for better visualization
axis equal;
