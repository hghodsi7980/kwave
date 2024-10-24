% Parameters
numPoints = 5000;
a = round(max(20, abs(normrnd(75, 25))));
b = round(max(20, abs(normrnd(75, 25))));
c = round(max(20, abs(normrnd(75, 25))));

% Generate random points on ellipsoid
theta = 2 * pi * rand(numPoints, 1);
phi = pi * rand(numPoints, 1);
x_ellipsoid = a * sin(phi) .* cos(theta);
y_ellipsoid = b * sin(phi) .* sin(theta);
z_ellipsoid = c * cos(phi);

% Generate ellipsoid points with random displacements
x_displacement = a * rand(numPoints, 1);
y_displacement = b * rand(numPoints, 1);
z_displacement = c * rand(numPoints, 1);
x = x_ellipsoid + x_displacement;
y = y_ellipsoid + y_displacement;
z = z_ellipsoid + z_displacement;

% Create Delaunay triangulation
DT = delaunayTriangulation(x, y, z);

% Plot the surface of the ellipsoid
view(125, 25);
axis equal;
hold on;

% Generate random points for clusters
min_points = 4000;
num_clusters = max(5, randi(50));
number_of_overall_points = round(max(min_points, abs(normrnd(10000, 2000))));
cluster_centers = DT.Points;
num_points_to_exclude = size(DT.Points, 1) - num_clusters;
indices_to_exclude = randperm(size(cluster_centers, 1), num_points_to_exclude);
cluster_centers(indices_to_exclude, :) = [];

% Define bandwidth for Gaussian kernel
bandwidth = 600 / num_clusters;
coordinates_to_check = [];

% Generate points using KDE around each cluster center
for i = 1:num_clusters
    cluster_points = mvnrnd(cluster_centers(i, :), bandwidth^2 * eye(3), round(number_of_overall_points / num_clusters));
    coordinates_to_check = [coordinates_to_check; cluster_points];
end

% Check if points are inside or on the surface of the ellipsoid
in_or_on_ellipsoid = ~isnan(pointLocation(DT, coordinates_to_check));

% Select inside points
inside_points = coordinates_to_check(in_or_on_ellipsoid, :);
Points = unique(inside_points, 'rows');

% Plot the inside points
scatter3(Points(:, 1), Points(:, 2), Points(:, 3), '.', 'b');

% Initialize an array to store the number of connections for each point
num_connections = zeros(size(Points, 1), 1);
max_connections = 4;  % Each point can have up to 4 connections

% Initialize variable to store all points including intermediate ones
all_points = Points;
adjacency_list = {}; % Initialize adjacency list to keep track of connections

% Function to calculate the angle between two vectors
calculate_angle = @(v1, v2) acos(dot(v1, v2) / (norm(v1) * norm(v2)));

% Find the 3 or 4 nearest neighbors for each point
for i = 1:size(Points, 1)
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
            
            % Control angles between connections
            if connected_neighbors > 1
                valid_angle = true;
                for k = 1:connected_neighbors
                    existing_neighbor_idx = neighbors_idx(k);
                    existing_bond = Points(existing_neighbor_idx, :) - point;
                    new_bond = neighbor_point - point;
                    
                    % Calculate angle between existing bond and new bond
                    angle_between = calculate_angle(existing_bond, new_bond);
                    
                    % If angle is less than 30 degrees, skip this connection
                    if angle_between < pi / 6  % 30 degrees in radians
                        valid_angle = false;
                        break;
                    end
                end
                if ~valid_angle
                    continue;  % Skip this connection
                end
            end
            
            % If bond length exceeds 15, add intermediate points
            if bond_length > 15
                num_intermediate_points = ceil(bond_length / 15);  % Calculate how many points are needed
                step_vector = (neighbor_point - point) / num_intermediate_points;  % Vector for each step
                
                % Add intermediate points along the bond
                previous_point = point; % Store the previous point
                for k = 1:num_intermediate_points
                    new_point = point + k * step_vector;
                    all_points = [all_points; new_point];  % Append new intermediate point to all points
                    adjacency_list{end+1} = [previous_point; new_point]; % Add to adjacency list
                    previous_point = new_point; % Update previous point to new point
                end
                
                % Add final connection to the original neighbor
                adjacency_list{end+1} = [previous_point; neighbor_point]; % Add the final connection to the list
            else
                % Directly connect the points if bond length is within the limit
                adjacency_list{end+1} = [point; neighbor_point]; % Add connection to adjacency list
            end
            
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

% Now plot the connections from the adjacency list, while ensuring that only valid points are connected
for i = 1:length(adjacency_list)
    connection = adjacency_list{i};
    % Check if the length of the connection exceeds 15 units
    bond_length = norm(connection(1, :) - connection(2, :));
    if bond_length <= 15  % Plot only if the bond length is within the limit
        plot3([connection(1, 1) connection(2, 1)], [connection(1, 2) connection(2, 2)], [connection(1, 3) connection(2, 3)], 'Color', [1 0 0 0.4]);
    end
end

% Check for isolated points and ensure full connectivity
for i = 1:size(Points, 1)
    if num_connections(i) == 0
        % If any points are isolated, connect them to their nearest neighbor
        nearest_neighbor_idx = knnsearch(Points, Points(i, :), 'K', 2); % Find nearest neighbor
        nearest_neighbor_idx = nearest_neighbor_idx(2);  % Ignore self (first element)
        neighbor_point = Points(nearest_neighbor_idx, :);
        adjacency_list{end+1} = [Points(i, :); neighbor_point]; % Add the connection to the list
        plot3([Points(i, 1) neighbor_point(1)], [Points(i, 2) neighbor_point(2)], [Points(i, 3) neighbor_point(3)], 'Color', [1 0 0 0.4]);
    end
end

hold off;
