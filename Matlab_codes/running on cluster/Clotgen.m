function Clotgen(startIter, endIter)
    % If not provided, use default values for startIter and endIter
    if nargin < 2
        startIter = 1;
        endIter = 10;
    end

    for index = startIter:endIter
    % Set the random seed for reproducibility
    fixed_seed = 5000;  % Choose a fixed seed value
    rng(fixed_seed + index);  % Randomize the seed in each iteration based on index

    % Parameters
    numPoints = 5000; % Number of points on ellipsoid surface
    a = round(max(20, abs(normrnd(75, 25))));
    b = round(max(20, abs(normrnd(75, 25))));
    c = round(max(20, abs(normrnd(75, 25))));

    % Random displacement factors
    displacementFactor_x = 1 + 0.5 * rand();
    displacementFactor_y = 1 + 0.5 * rand();
    displacementFactor_z = 1 + 0.5 * rand();

    % Generate random points on ellipsoid surface
    theta = 2 * pi * rand(numPoints, 1);
    phi = pi * rand(numPoints, 1);
    x_ellipsoid = a * sin(phi) .* cos(theta);
    y_ellipsoid = b * sin(phi) .* sin(theta);
    z_ellipsoid = c * cos(phi);

    % Add random displacement
    x_displacement = displacementFactor_x * a * rand(numPoints, 1);
    y_displacement = displacementFactor_y * b * rand(numPoints, 1);
    z_displacement = displacementFactor_z * c * rand(numPoints, 1);
    x = x_ellipsoid + x_displacement;
    y = y_ellipsoid + y_displacement;
    z = z_ellipsoid + z_displacement;

    % Create Delaunay triangulation
    DT = delaunayTriangulation(x, y, z);

    % Calculate convex hull
    [C, volume] = convexHull(DT);

    % Generate random points for clusters
    min_points = 500;
    num_clusters = max(5, randi(50));
    number_of_overall_points = round(max(min_points, abs(normrnd(2000, 1000))));
    cluster_centers = DT.Points;
    num_points_to_exclude = size(DT.Points, 1) - num_clusters;
    indices_to_exclude = randperm(size(cluster_centers, 1), num_points_to_exclude);
    cluster_centers(indices_to_exclude, :) = [];

    % Define bandwidth for Gaussian kernel
    bandwidth = 1000 / num_clusters;
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

    % Create Delaunay triangulation for inside points
    DT_inside = delaunayTriangulation(Points);
    [V, R] = voronoiDiagram(DT_inside);

    % Remove rows with Inf values from V
    finite_indices = all(isfinite(V), 2);
    V_finite = V(finite_indices, :);

    % Mapping from old indices to new indices
    old_to_new_index = zeros(size(V, 1), 1);
    old_to_new_index(finite_indices) = 1:sum(finite_indices);

    % Determine points inside the convex hull
    hull_tetrahedrons = delaunayTriangulation(DT.Points(C(:), :));
    inside_hull = ~isnan(tsearchn(hull_tetrahedrons.Points, hull_tetrahedrons.ConnectivityList, V_finite));

    % Update vertices inside the convex hull
    V_finite_inside = V_finite(inside_hull, :);
    old_to_new_inside_index = zeros(size(V_finite, 1), 1);
    old_to_new_inside_index(inside_hull) = 1:sum(inside_hull);

    % Create edges for Voronoi diagram
    edges = [];
    for i = 1:length(R)
        region = R{i};
        finite_region = old_to_new_index(region(region <= size(V, 1) & region > 0));
        finite_region = finite_region(finite_region > 0);
        finite_region = old_to_new_inside_index(finite_region);
        finite_region = finite_region(finite_region > 0);
        if numel(finite_region) > 1
            region_edges = [finite_region(:), circshift(finite_region(:), -1)];
            edges = [edges; region_edges];
        end
    end
    edges = unique(sort(edges, 2), 'rows');
    Points = V_finite_inside;

    % Calculate edge lengths
    edge_lengths = sqrt(sum((V_finite_inside(edges(:, 1), :) - V_finite_inside(edges(:, 2), :)).^2, 2));

    % Set distance threshold for edge removal
    distance_threshold = 1;
    remove_probabilities = 1 - exp(-edge_lengths / distance_threshold);

    % Initialize graph and calculate centrality
    G = graph(edges(:, 1), edges(:, 2));
    CC = conncomp(G);
    degree_cent = centrality(G, 'degree');
    mean_rank = mean(degree_cent);
    mean_rank_best = 3.5;
    std_rank = std(degree_cent);
    num_edges_init = size(edges, 1);

    % Set angle thresholds
    low_angle_threshold = pi / 6;
    high_angle_threshold = 3 * pi / 2;

    % Outer loop to ensure the mean rank is below the threshold
    while true
        G = graph(edges(:, 1), edges(:, 2));
        degree_cent = centrality(G, 'degree');
        mean_rank = mean(degree_cent);

        for node = 1:size(Points, 1)
            force_remove = 0;
            connected_edges = edges(any(edges == node, 2), :);
            if size(connected_edges, 1) > 4
                force_remove = 1;
            end

            for i = 1:size(connected_edges, 1)
                edge = connected_edges(i, :);
                ind = find(ismember(edges, edge, 'rows'), 1);

                if force_remove && remove_probabilities(ind) > 0.5
                    remove_probabilities(ind) = 1;
                end

                if rand() < remove_probabilities(ind)
                    % Identify affected nodes
                    affected_nodes = edge;
                    new_edges = edges([1:ind-1, ind+1:end], :);
                    G_i = graph(new_edges(:, 1), new_edges(:, 2));
                    CC_i = conncomp(G_i);
                    degree_cent_i = centrality(G_i, 'degree');
                    mean_rank_i = mean(degree_cent_i);
                    std_rank_i = std(degree_cent_i);

                    % Calculate angles for affected nodes
                    angles = calculate_angles_for_nodes(Points, connected_edges);

                    % Check constraints
                    if max(CC) == max(CC_i) && min(degree_cent_i) > 2 && ...
                            (abs(mean_rank_i - mean_rank_best) <= abs(mean_rank - mean_rank_best) || std_rank_i < std_rank)
                        if any(angles < low_angle_threshold) || any(angles > high_angle_threshold)
                            node_edges = edges(any(edges == node, 2), :);
                            node_edge_inds = find(any(ismember(edges, node_edges), 2));
                            max_prob = max(remove_probabilities(node_edge_inds));
                            remove_probabilities(node_edge_inds) = remove_probabilities(node_edge_inds) * 2;
                            if force_remove
                                remove_probabilities(node_edge_inds) = remove_probabilities(node_edge_inds) * 5;
                            end
                        else
                            % Update the graph if constraints are met
                            edges = new_edges;
                            G = graph(edges(:, 1), edges(:, 2));
                            CC = conncomp(G);
                            degree_cent = centrality(G, 'degree');
                            mean_rank = mean(degree_cent);
                            std_rank = std(degree_cent);
                            % disp(100 * mean_rank_best / mean_rank)
                        end
                    end
                end
            end
        end

        % Check if the mean rank is below the threshold
        if mean_rank < mean_rank_best
            break;
        end

        % Adjust probabilities
        remove_probabilities = remove_probabilities * 2;
        remove_probabilities(remove_probabilities > 1) = 1;
        low_angle_threshold = low_angle_threshold * 0.5;
        high_angle_threshold = high_angle_threshold * 1.5;
    end

    % Update edges to make outer edges realistic
    [~, edges] = find_outer_points_and_update_edges(V_finite_inside, edges);

    % Initialize coordinates and material indices
    all_coordinates = [];
    line_coordinates = [];
    % hold on;
    % axis equal;

    % Exclude some points randomly
    exclude_percentage_platelet = rand();
    num_points_to_exclude = round(exclude_percentage_platelet * size(Points, 1));
    indices_to_exclude = randperm(size(Points, 1), num_points_to_exclude);
    Points(indices_to_exclude, :) = [];

    % Initialize platelet parameters
    platelet_radius = 1;
    for i = 1:size(Points, 1)
        % Generate and plot platelet spheres
        [sx, sy, sz] = sphere;
        sx = sx * platelet_radius + Points(i, 1);
        sy = sy * platelet_radius + Points(i, 2);
        sz = sz * platelet_radius + Points(i, 3);
        % h = surf(sx, sy, sz, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        % set(h, 'FaceColor', [0 0.8 0.13]);
        % pause(1e-12);

        % Append surface sphere coordinates
        sphere_surface_coordinates = [sx(:), sy(:), sz(:), ones(size(sx(:))) * 3];
        all_coordinates = [all_coordinates; sphere_surface_coordinates];

        % Generate inside sphere coordinates
        num_inside_points = 250;
        sphere_center = Points(i, :);
        sphere_radius = platelet_radius;
        rand_points = rand(num_inside_points, 3) - 0.5;
        rand_points = rand_points ./ sqrt(sum(rand_points.^2, 2));
        inside_sphere_coordinates = sphere_center + sphere_radius * rand_points;
        sphere_inside_coordinates = [inside_sphere_coordinates, ones(size(inside_sphere_coordinates, 1), 1) * 3];
        all_coordinates = [all_coordinates; sphere_inside_coordinates];
    end

    for i = 1:size(edges, 1)
        start_point = V_finite_inside(edges(i, 1), :);
        end_point = V_finite_inside(edges(i, 2), :);
        interp_pts = [];
        for t = linspace(0, 1, 100)
            interp_pt = (1 - t) * start_point + t * end_point;
            interp_pts = [interp_pts; interp_pt];
        end
        line_coordinates = [line_coordinates; interp_pts];
        % h = plot3(interp_pts(:, 1), interp_pts(:, 2), interp_pts(:, 3), 'm', 'LineWidth', 1);
        % set(h, 'Color', [0 0.4 0.13 0.5]);
        % pause(1e-12);
    end

    material_index_line = ones(size(line_coordinates, 1), 1) * 2;
    line_coordinates = [line_coordinates, material_index_line];
    all_coordinates = [all_coordinates; line_coordinates];

    % RBC filling factor
    RBC_filling_factor = 0.7+0.3*rand();
    RBC_diameter = 8;
    for i = 1:size(inside_points, 1)
        pause(1e-12);
        if rand() < RBC_filling_factor
            center_point_t = inside_points(i, :);
            distances = zeros(round(size(line_coordinates, 1) / 10), 1);
            for p = 1:10:size(line_coordinates, 1)
                distances(round(p / 10) + 1) = sqrt(sum((line_coordinates(p, 1:3) - center_point_t).^2));
            end
            distances_sorted = sort(distances);
            if distances_sorted(1) > RBC_diameter && distances_sorted(1) < 1.3 * RBC_diameter
                scaling_factor_RBC = 1;
                cell_rotations = rand(1, 3) * pi;
                d = 8 / scaling_factor_RBC;
                br = 1 / scaling_factor_RBC;
                h = 2.12 / scaling_factor_RBC;
                P = -(d^2 / 2) + (h^2 / 2) * ((d^2 / br^2) - 1) - h^2 / 2 * ((d^2 / br^2) - 1) * sqrt(1 - (br^2 / h^2));
                Q = P * (d^2 / br^2) + (br^2 / 4) * (d^4 / br^4 - 1);
                R = -P * (d^2 / 4) - d^4 / 16;
                [x_rbc, y_rbc, z_rbc] = meshgrid(-10 + center_point_t(1):0.5:10 + center_point_t(1), -10 + center_point_t(2):0.5:10 + center_point_t(2), -10 + center_point_t(3):0.5:10 + center_point_t(3));
                x_rbc_rot = x_rbc - center_point_t(1);
                y_rbc_rot = y_rbc - center_point_t(2);
                z_rbc_rot = z_rbc - center_point_t(3);
                x_rbc_temp = x_rbc_rot;
                x_rbc_rot = x_rbc_temp * cos(cell_rotations(1)) - z_rbc_rot * sin(cell_rotations(1));
                z_rbc_rot = x_rbc_temp * sin(cell_rotations(1)) + z_rbc_rot * cos(cell_rotations(1));
                y_rbc_temp = y_rbc_rot;
                y_rbc_rot = y_rbc_temp * cos(cell_rotations(2)) + z_rbc_rot * sin(cell_rotations(2));
                z_rbc_rot = -y_rbc_temp * sin(cell_rotations(2)) + z_rbc_rot * cos(cell_rotations(2));
                x_rbc_temp = x_rbc_rot;
                x_rbc_rot = x_rbc_temp * cos(cell_rotations(3)) - y_rbc_rot * sin(cell_rotations(3));
                y_rbc_rot = x_rbc_temp * sin(cell_rotations(3)) + y_rbc_rot * cos(cell_rotations(3));
                eq = ((x_rbc_rot).^2 + (y_rbc_rot).^2 + (z_rbc_rot).^2).^2 + P * ((x_rbc_rot).^2 + (y_rbc_rot).^2) + Q * (z_rbc_rot).^2 + R;
                v_rbc = zeros(size(x_rbc));
                v_rbc(eq <= 0) = 1;
                rbc_x = x_rbc(v_rbc == 1);
                rbc_y = y_rbc(v_rbc == 1);
                rbc_z = z_rbc(v_rbc == 1);
                rbc_coordinates = [rbc_x(:), rbc_y(:), rbc_z(:)];
                % DTr = delaunayTriangulation(rbc_x, rbc_y, rbc_z);
                % [Cr, vr] = convexHull(DTr);
                % trisurf(Cr, DTr.Points(:, 1), DTr.Points(:, 2), DTr.Points(:, 3), 'FaceColor', 'red', 'EdgeColor', 'none');
                material_index_rbc = ones(size(rbc_coordinates, 1), 1) * 1;
                all_coordinates = [all_coordinates; rbc_coordinates, material_index_rbc];
            elseif distances_sorted(1) > 1.3 * RBC_diameter
                num_RBC_in_region = round(10 * RBC_filling_factor * distances_sorted(1) / RBC_diameter);
                for j = 1:num_RBC_in_region
                    intersect = true;
                    while intersect
                        center_point_t_1 = center_point_t + ((-1).^randi(2,1,3)) .* rand(1,3) * distances_sorted(1);
                        scaling_factor_RBC = 1;
                        cell_rotations = rand(1, 3) * pi;
                        d = 8 / scaling_factor_RBC;
                        br = 1 / scaling_factor_RBC;
                        h = 2.12 / scaling_factor_RBC;
                        P = -(d^2 / 2) + (h^2 / 2) * ((d^2 / br^2) - 1) - h^2 / 2 * ((d^2 / br^2) - 1) * sqrt(1 - (br^2 / h^2));
                        Q = P * (d^2 / br^2) + (br^2 / 4) * (d^4 / br^4 - 1);
                        R = -P * (d^2 / 4) - d^4 / 16;
                        [x_rbc, y_rbc, z_rbc] = meshgrid(-10 + center_point_t_1(1):0.5:10 + center_point_t_1(1), -10 + center_point_t_1(2):0.5:10 + center_point_t_1(2), -10 + center_point_t_1(3):0.5:10 + center_point_t_1(3));
                        x_rbc_rot = x_rbc - center_point_t_1(1);
                        y_rbc_rot = y_rbc - center_point_t_1(2);
                        z_rbc_rot = z_rbc - center_point_t_1(3);
                        x_rbc_temp = x_rbc_rot;
                        x_rbc_rot = x_rbc_temp * cos(cell_rotations(1)) - z_rbc_rot * sin(cell_rotations(1));
                        z_rbc_rot = x_rbc_temp * sin(cell_rotations(1)) + z_rbc_rot * cos(cell_rotations(1));
                        y_rbc_temp = y_rbc_rot;
                        y_rbc_rot = y_rbc_temp * cos(cell_rotations(2)) + z_rbc_rot * sin(cell_rotations(2));
                        z_rbc_rot = -y_rbc_temp * sin(cell_rotations(2)) + z_rbc_rot * cos(cell_rotations(2));
                        x_rbc_temp = x_rbc_rot;
                        x_rbc_rot = x_rbc_temp * cos(cell_rotations(3)) - y_rbc_rot * sin(cell_rotations(3));
                        y_rbc_rot = x_rbc_temp * sin(cell_rotations(3)) + y_rbc_rot * cos(cell_rotations(3));
                        eq = ((x_rbc_rot).^2 + (y_rbc_rot).^2 + (z_rbc_rot).^2).^2 + P * ((x_rbc_rot).^2 + (y_rbc_rot).^2) + Q * (z_rbc_rot).^2 + R;
                        v_rbc = zeros(size(x_rbc));
                        v_rbc(eq <= 0) = 1;
                        rbc_x = x_rbc(v_rbc == 1);
                        rbc_y = y_rbc(v_rbc == 1);
                        rbc_z = z_rbc(v_rbc == 1);
                        rbc_coordinates = [rbc_x(:), rbc_y(:), rbc_z(:)];
                        min_d = 1.0;
                        squared_distances_to_center = sum((line_coordinates(:, 1:3) - center_point_t_1).^2, 2);
                        inside_sphere = line_coordinates(squared_distances_to_center <= (RBC_diameter / 2)^2, 1:3);
                        distances_to_Fibrin = distances_sorted(1);
                        if ~isempty(inside_sphere)
                            distances_to_Fibrin = 2 * min_d;
                        else
                            for k = 1:size(rbc_coordinates, 1)
                                for l = 1:size(inside_sphere, 1)
                                    distances_to_Fibrin = sqrt(sum((rbc_coordinates(k, 1:3) - inside_sphere(l, :)).^2));
                                    if distances_to_Fibrin <= min_d
                                        break;
                                    end
                                end
                            end
                        end
                        if distances_to_Fibrin > min_d
                            % DTr = delaunayTriangulation(rbc_x, rbc_y, rbc_z);
                            % [Cr, vr] = convexHull(DTr);
                            % trisurf(Cr, DTr.Points(:, 1), DTr.Points(:, 2), DTr.Points(:, 3), 'FaceColor', 'red', 'EdgeColor', 'none');
                            material_index_rbc = ones(size(rbc_coordinates, 1), 1) * 1;
                            all_coordinates = [all_coordinates; rbc_coordinates, material_index_rbc];
                            intersect = false;
                        else
                            intersect = true;
                            scaling_factor_RBC = scaling_factor_RBC * sqrt(0.9);
                        end
                    end
                end
            elseif distances_sorted(1) <= RBC_diameter && distances_sorted(1) > RBC_diameter / 4
                scaling_factor_RBC = sqrt(RBC_diameter / distances_sorted(1));
                cell_rotations = rand(1, 3) * pi;
                d = 8 / scaling_factor_RBC;
                br = 1 / scaling_factor_RBC;
                h = 2.12 / scaling_factor_RBC;
                P = -(d^2 / 2) + (h^2 / 2) * ((d^2 / br^2) - 1) - h^2 / 2 * ((d^2 / br^2) - 1) * sqrt(1 - (br^2 / h^2));
                Q = P * (d^2 / br^2) + (br^2 / 4) * (d^4 / br^4 - 1);
                R = -P * (d^2 / 4) - d^4 / 16;
                [x_rbc, y_rbc, z_rbc] = meshgrid(-10 + center_point_t(1):0.5:10 + center_point_t(1), -10 + center_point_t(2):0.5:10 + center_point_t(2), -10 + center_point_t(3):0.5:10 + center_point_t(3));
                x_rbc_rot = x_rbc - center_point_t(1);
                y_rbc_rot = y_rbc - center_point_t(2);
                z_rbc_rot = z_rbc - center_point_t(3);
                x_rbc_temp = x_rbc_rot;
                x_rbc_rot = x_rbc_temp * cos(cell_rotations(1)) - z_rbc_rot * sin(cell_rotations(1));
                z_rbc_rot = x_rbc_temp * sin(cell_rotations(1)) + z_rbc_rot * cos(cell_rotations(1));
                y_rbc_temp = y_rbc_rot;
                y_rbc_rot = y_rbc_temp * cos(cell_rotations(2)) + z_rbc_rot * sin(cell_rotations(2));
                z_rbc_rot = -y_rbc_temp * sin(cell_rotations(2)) + z_rbc_rot * cos(cell_rotations(2));
                x_rbc_temp = x_rbc_rot;
                x_rbc_rot = x_rbc_temp * cos(cell_rotations(3)) - y_rbc_rot * sin(cell_rotations(3));
                y_rbc_rot = x_rbc_temp * sin(cell_rotations(3)) + y_rbc_rot * cos(cell_rotations(3));
                eq = ((x_rbc_rot).^2 + (y_rbc_rot).^2 + (z_rbc_rot).^2).^2 + P * ((x_rbc_rot).^2 + (y_rbc_rot).^2) + Q * (z_rbc_rot).^2 + R;
                v_rbc = zeros(size(x_rbc));
                v_rbc(eq <= 0) = 1;
                rbc_x = x_rbc(v_rbc == 1);
                rbc_y = y_rbc(v_rbc == 1);
                rbc_z = z_rbc(v_rbc == 1);
                rbc_coordinates = [rbc_x(:), rbc_y(:), rbc_z(:)];
                % DTr = delaunayTriangulation(rbc_x, rbc_y, rbc_z);
                % [Cr, vr] = convexHull(DTr);
                % trisurf(Cr, DTr.Points(:, 1), DTr.Points(:, 2), DTr.Points(:, 3), 'FaceColor', 'red', 'EdgeColor', 'none');
                material_index_rbc = ones(size(rbc_coordinates, 1), 1) * 1;
                all_coordinates = [all_coordinates; rbc_coordinates, material_index_rbc];
            end
        end
    end

    % Sort the coordinates matrix based on first x, then y, and finally z
    all_coordinates_sorted = sortrows(single(all_coordinates), [4, 1, 2, 3]);
    coordinates = all_coordinates_sorted(:, 1:3);
    materials = all_coordinates_sorted(:, 4);

    % Define the size of the mesh space
    distance_to_boundary = 50;
    mesh_size_x = round(max(coordinates(:, 1)) - min(coordinates(:, 1))) + 2 * distance_to_boundary;
    mesh_size_y = round(max(coordinates(:, 2)) - min(coordinates(:, 2))) + 2 * distance_to_boundary;
    mesh_size_z = round(max(coordinates(:, 3)) - min(coordinates(:, 3))) + 2 * distance_to_boundary;

    % Create an empty mesh space
    mesh_space = zeros(mesh_size_x, mesh_size_y, mesh_size_z, 'single');

    % Translate coordinates to positive quadrant
    coordinates = coordinates - min(coordinates) + 1;

    % Scale coordinates to fit within the mesh size
    coordinates_scaled = zeros(size(coordinates, 1), 3); % Initialize here
    scaling_factor_x = 1;
    scaling_factor_y = 1;
    scaling_factor_z = 1;
    coordinates_scaled(:, 1) = coordinates(:, 1) * scaling_factor_x + distance_to_boundary;
    coordinates_scaled(:, 2) = coordinates(:, 2) * scaling_factor_y + distance_to_boundary;
    coordinates_scaled(:, 3) = coordinates(:, 3) * scaling_factor_z + distance_to_boundary;

    % Ensure all coordinates are integers
    coordinates_scaled = round(coordinates_scaled);

    % Assign material indices to the mesh space
    for i = 1:size(coordinates_scaled, 1)
        x = coordinates_scaled(i, 1);
        y = coordinates_scaled(i, 2);
        z = coordinates_scaled(i, 3);
        if x >= 1 && x <= mesh_size_x && y >= 1 && y <= mesh_size_y && z >= 1 && z <= mesh_size_z
            mesh_space(x, y, z) = materials(i);
        end
    end

    % Save the mesh space to a file
    save(sprintf('mesh_space_%d.mat', index), 'mesh_space', 'volume', '-v7.3');
    end
end

function [outer_points, all_edges] = find_outer_points_and_update_edges(points, edges)
    % Find outer points using convex hull
    K = convhull(points(:, 1), points(:, 2), points(:, 3));
    outer_point_indices = unique(K(:));
    outer_points = points(outer_point_indices, :);

    % Find edges connected to outer points
    outer_edges = [];
    for i = 1:size(edges, 1)
        if ismember(edges(i, 1), outer_point_indices) || ismember(edges(i, 2), outer_point_indices)
            outer_edges = [outer_edges; edges(i, :)];
        end
    end

    % Delete all but one edge connected to each outer point randomly
    updated_edges = outer_edges;
    outer_point_used = zeros(size(outer_point_indices));
    for i = 1:numel(outer_point_indices)
        idx = find(outer_edges(:, 1) == outer_point_indices(i) | outer_edges(:, 2) == outer_point_indices(i));
        if numel(idx) > 1
            idx = idx(randperm(numel(idx)));
            if numel(idx) > 1
                updated_edges(idx(2:end), :) = ones(length(idx(2:end)), 2);
            end
        end
    end

    % Construct all_edges matrix including updated_edges and other edges
    all_edges = edges;
    [~, ia, ~] = intersect(all_edges, updated_edges, 'rows');
    all_edges(ia, :) = [];
end

function angles = calculate_angles_for_nodes(Points, connected_edges)
    % Calculate vectors based on connected edges
    vectors = Points(connected_edges(:, 1), :) - Points(connected_edges(:, 2), :);
    vectors = vectors ./ vecnorm(vectors, 2, 2);

    % Initialize angles array
    angles = zeros(size(vectors, 1), 1);

    % Calculate angles between consecutive vectors
    for i = 1:size(vectors, 1)
        vec1 = vectors(i, :);
        if i < size(vectors, 1)
            vec2 = vectors(i + 1, :);
        else
            vec2 = vectors(1, :);
        end
        dot_product = dot(vec1, vec2);
        angles(i) = acos(dot_product);
    end
end
