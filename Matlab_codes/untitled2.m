% Parameters
numPoints = 5000;
a = round(max(20, abs(normrnd(75, 25))));
b = round(max(20, abs(normrnd(75, 25))));
c = round(max(20, abs(normrnd(75, 25))));

% Generate random points uniformly on ellipsoid surface
theta = 2 * pi * rand(numPoints, 1);
phi = pi * rand(numPoints, 1);
x_ellipsoid = a * sin(phi) .* cos(theta);
y_ellipsoid = b * sin(phi) .* sin(theta);
z_ellipsoid = c * cos(phi);

% Add uniform random displacements
x_displacement = a * (rand(numPoints, 1) - 0.5) * 2;
y_displacement = b * (rand(numPoints, 1) - 0.5) * 2;
z_displacement = c * (rand(numPoints, 1) - 0.5) * 2;
x = x_ellipsoid + x_displacement;
y = y_ellipsoid + y_displacement;
z = z_ellipsoid + z_displacement;

% Store points
Points = [x, y, z];

% Plot the initial points
figure;
scatter3(Points(:, 1), Points(:, 2), Points(:, 3), 'b.');
hold on;
axis equal;

% Initialize an array to store the number of connections for each point
num_connections = zeros(numPoints, 1);
max_connections = 4;  % Each point can have up to 4 connections

% Define bond length limit
max_bond_length = 15;

% Find nearest neighbors and apply constraints
for i = 1:numPoints
    point = Points(i, :);

    % Find nearest neighbors excluding the point itself
    neighbors_idx = knnsearch(Points, point, 'K', max_connections + 1); % K+1 to include the point itself
    neighbors_idx(1) = []; % Remove the point itself from neighbors

    connected_neighbors = 0;  % Count of connected neighbors for current point

    for j = 1:length(neighbors_idx)
        neighbor_idx = neighbors_idx(j);

        % Check if the point and the neighbor already have enough connections
        if num_connections(i) < max_connections && num_connections(neighbor_idx) < max_connections
            % Check bond length
            neighbor_point = Points(neighbor_idx, :);
            bond_length = norm(point - neighbor_point);

            % If bond length exceeds the limit, move the neighbor closer
            if bond_length > max_bond_length
                direction = (neighbor_point - point) / bond_length; % Unit vector in bond direction
                new_position = point + max_bond_length * direction;  % Move neighbor to within the bond length limit
                Points(neighbor_idx, :) = new_position;  % Update neighbor's position
                neighbor_point = new_position;
                bond_length = max_bond_length;  % Adjust bond length to the maximum allowed
            end

            % Plot the connection as a red line
            plot3([point(1), neighbor_point(1)], [point(2), neighbor_point(2)], [point(3), neighbor_point(3)], 'r-');

            % Update connection counts for both the point and the neighbor
            num_connections(i) = num_connections(i) + 1;
            num_connections(neighbor_idx) = num_connections(neighbor_idx) + 1;
            connected_neighbors = connected_neighbors + 1;
        end

        % Stop connecting if the current point reaches the max number of connections
        if num_connections(i) >= max_connections
            break;
        end
    end
end

% Check for isolated points and ensure full connectivity
for i = 1:numPoints
    if num_connections(i) == 0
        % If any points are isolated, connect them to their nearest neighbor
        nearest_neighbor_idx = knnsearch(Points, Points(i, :), 'K', 2); % Find nearest neighbor
        nearest_neighbor_idx = nearest_neighbor_idx(2);  % Ignore self (first element)
        neighbor_point = Points(nearest_neighbor_idx, :);
        plot3([Points(i, 1), neighbor_point(1)], [Points(i, 2), neighbor_point(2)], [Points(i, 3), neighbor_point(3)], 'g-');
    end
end

hold off;
