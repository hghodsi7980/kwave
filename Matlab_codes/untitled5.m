% Example 3D coordinates
coordinates = [
    1, 2, 3;
    4, 5, 6;
    7, 8, 9;
    % Add more coordinates as needed
];

% Calculate the pairwise distances
num_points = size(coordinates, 1);
distances = zeros(num_points, num_points);

for i = 1:num_points
    for j = i+1:num_points
        distances(i, j) = sqrt(sum((coordinates(i,:) - coordinates(j,:)).^2));
        distances(j, i) = distances(i, j); % Symmetric matrix
    end
end

% Find immediate neighbors
mean_distances = zeros(1, num_points);
for i = 1:num_points
    neighbor_distances = distances(i, :);
    neighbor_distances(i) = NaN; % Ignore distance to itself
    neighbor_distances(isnan(neighbor_distances)) = []; % Remove NaN values
    mean_distances(i) = mean(neighbor_distances(1:min(2, numel(neighbor_distances)))); % Considering only two closest neighbors
end

% Calculate mean distance between adjacent points
mean_distance = mean(mean_distances(~isnan(mean_distances)));

disp(['Mean distance between adjacent points: ', num2str(mean_distance)]);
