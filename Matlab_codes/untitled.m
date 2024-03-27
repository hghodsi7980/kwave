


% Ensure connectivity using connected components
G = graph(remaining_edges(:,1), remaining_edges(:,2));

% Find connected components
CC = conncomp(G);

% Check if all points are already connected
if max(CC) == 1
    fully_connected_edges = remaining_edges;
else
    % Connect the components
    while max(CC) > 1
        % Find points in different components
        component_indices = unique(CC);
        num_components = length(component_indices);
        if num_components > 1
            % Select one point from each component
            point_indices = zeros(num_components, 1);
            for i = 1:num_components
                point_indices(i) = find(CC == component_indices(i), 1);
            end

            % Find the closest pair of points from different components
            min_distance = inf;
            closest_pair = [];
            for i = 1:num_components-1
                for j = i+1:num_components
                    distance_ij = norm(points(point_indices(i), :) - points(point_indices(j), :));
                    if distance_ij < min_distance
                        min_distance = distance_ij;
                        closest_pair = [point_indices(i), point_indices(j)];
                    end
                end
            end

            % Add edge between the closest pair
            remaining_edges = [remaining_edges; closest_pair];

            % Update connected components
            G = graph(remaining_edges(:,1), remaining_edges(:,2));
            CC = conncomp(G);
        end
    end

    % Extract edges of the fully connected graph
    fully_connected_edges = remaining_edges;
end

% Plot the fully connected edges using plot3
figure;
hold on;
for i = 1:size(fully_connected_edges, 1)
    start_point = points(fully_connected_edges(i, 1), :);
    end_point = points(fully_connected_edges(i, 2), :);
    plot3([start_point(1), end_point(1)], [start_point(2), end_point(2)], [start_point(3), end_point(3)], 'b');
end
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Fully Connected Edges');
grid on;
hold off;
