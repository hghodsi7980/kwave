for index = 1:1
    clearvars('-except', 'index');
    close all
    hold on
    grid on
    grid minor
    axis equal
    % Parameters
    numPoints = 5000; % Adjust as needed

    % Define clot size um
    a = round(max(200, normrnd(250,75))); % Semi-axis length along x-axis
    b = round(max(200, normrnd(250,75))); % Semi-axis length along y-axis
    c = round(max(200, normrnd(250,75))); % Semi-axis length along z-axis

    displacementFactor = 1; % Control parameter for randomness

    % Generate random points on the surface of an ellipsoid
    theta = 2*pi*rand(numPoints, 1);
    phi = pi*rand(numPoints, 1);
    x_ellipsoid = a * sin(phi) .* cos(theta);
    y_ellipsoid = b * sin(phi) .* sin(theta);
    z_ellipsoid = c * cos(phi);

    % Add random displacement
    x_displacement = displacementFactor *a*(rand(numPoints, 1));
    y_displacement = displacementFactor *b*(rand(numPoints, 1));
    z_displacement = displacementFactor *c*(rand(numPoints, 1));
    x = x_ellipsoid + x_displacement;
    y = y_ellipsoid + y_displacement;
    z = z_ellipsoid + z_displacement;

    % Create Delaunay triangulation
    DT = delaunayTriangulation(x, y, z);

    % Calculate the convex hull
    [C,v] = convexHull(DT);

    % % trisurf(C,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3), ...
    % % 'EdgeColor','black','FaceColor','none')

    % minumum number of internal points  
    min_points = 1000; 

    % Set of coordinates to check
    number_of_overall_points = round(max(max([a,b,c])*4*abs((3+randn())),min_points));
    concentration_factor = 0.7;
    % Calculate the center point of the surface
    center_point = [mean(DT.Points(C, 1)), mean(DT.Points(C, 2)), mean(DT.Points(C, 3))];
    x_center = center_point(1);
    y_center = center_point(2);
    z_center = center_point(3);
    % The points on the surface
    Points = DT.Points(C,:);

    x_gaussian = normrnd(x_center, concentration_factor*a, number_of_overall_points, 1); % Adjust standard deviation as needed
    y_gaussian = normrnd(y_center, concentration_factor*b, number_of_overall_points, 1);
    z_gaussian = normrnd(z_center, concentration_factor*c, number_of_overall_points, 1);
    coordinates_to_check = [x_gaussian,y_gaussian,z_gaussian];
    % Check if points are inside or on the surface
    in_or_on_ellipsoid = ~isnan(pointLocation(DT, coordinates_to_check));

    % Choose the inside points
    inside_points = coordinates_to_check(in_or_on_ellipsoid, :);
    % Add inside points to the surface points
    Points = [Points;inside_points];
    Points = unique(Points, 'rows');
    % scatter3(Points(:,1),Points(:,2),Points(:,3),'filled','g')
    % Create Delaunay triangulation
    DT_inside = delaunayTriangulation(Points);

    % Get edge list
    edges = edges(DT_inside);

    % Calculate edge lengths
    edge_lengths = sqrt(sum((DT_inside.Points(edges(:, 1), :) - DT_inside.Points(edges(:, 2), :)).^2, 2));

    % Set distance threshold (adjust as needed)
    distance_threshold = 15; % You can adjust this threshold to control density

    % Calculate probability of edge removal based on distance
    remove_probabilities = 1 - exp(-edge_lengths / distance_threshold);
    G = graph(edges(:,1),edges(:,2));
    CC = conncomp(G);
    degree_cent = G.centrality("degree");
    mean_rank = mean(degree_cent);
    mean_rank_best = 3.5;
    std_rank = std(degree_cent);
    num_edges_init = size(edges,1);
    remaining_edges = edges;
    for i = 1:num_edges_init
        ind = find(remaining_edges== edges(i,:), 1 );
        if rand() < remove_probabilities(i)
            G_i = graph(remaining_edges([1:ind-1, ind+1:end],1),remaining_edges([1:ind-1, ind+1:end],2));
            CC_i = conncomp(G_i);
            degree_cent_i = G_i.centrality("degree");
            mean_rank_i = mean(degree_cent_i);
            std_rank_i = std(degree_cent_i);
            if max(CC) == max(CC_i) && (abs(mean_rank_i-mean_rank_best) <= abs(mean_rank-mean_rank_best) ...
                    || std_rank_i < std_rank ) && min(degree_cent_i) > 2
                remaining_edges = remaining_edges([1:ind-1, ind+1:end],:);
                G = graph(remaining_edges(:,1),remaining_edges(:,2));
                CC = conncomp(G);
                degree_cent = G.centrality("degree");
                mean_rank = mean(degree_cent);
                std_rank = std(degree_cent);
            end
        end
    end
    bond_count = size(remaining_edges,1);
    remaining_edges_test = remaining_edges;

    % Initialize matrix to store all coordinates and material indices
    all_coordinates = [];
    line_coordinates = []; % Initialize line coordinates
    rbc_coordinates = [];
    for i = 1:size(remaining_edges_test, 1)
        start_point = DT_inside.Points(remaining_edges_test(i, 1), :);
        end_point = DT_inside.Points(remaining_edges_test(i, 2), :);
        % Linearly interpolate points between pt1 and pt2
        interp_pts = [];
        for t = linspace(0, 1, 100)
            interp_pt = (1 - t) * start_point + t * end_point;
            interp_pts = [interp_pts; interp_pt];
        end
        % Append interpolated points to line_coordinates
        line_coordinates = [line_coordinates; interp_pts];
        plot3([start_point(1), end_point(1)], [start_point(2), end_point(2)], [start_point(3), end_point(3)], 'b');
        pause(1e-12);
    end
    % Generate material indices for each point
    material_index_line = ones(size(line_coordinates, 1), 1) * 2; % Material index for lines
    % Appending to line coordinates
    line_coordinates = [line_coordinates, material_index_line];
    % Append line coordinates to all_coordinates
    all_coordinates = [all_coordinates; line_coordinates];

    % Extract tetrahedra
    tetrahedra = DT_inside.ConnectivityList;

    % Initialize a list to store indices of tetrahedra with at least one remaining edge
    tetrahedra_with_remaining_edges = [];

    % Iterate through each tetrahedron
    for i = 1:size(tetrahedra, 1)
        tetrahedron_edges = [tetrahedra(i, [1 2]); tetrahedra(i, [1 3]); tetrahedra(i, [1 4]); ...
            tetrahedra(i, [2 3]); tetrahedra(i, [2 4]); tetrahedra(i, [3 4])];
        % Check if any edge of the tetrahedron is in the list of remaining edges
        if any(ismember(tetrahedron_edges, remaining_edges, 'rows'))
            tetrahedra_with_remaining_edges = [tetrahedra_with_remaining_edges; i];
        end
    end

    % Extract tetrahedra with at least one remaining edge
    tetrahedra_with_remaining_edges = tetrahedra(tetrahedra_with_remaining_edges, :);

    % Generate waitbar
    holder = waitbar(0, 'Progress RBC'); % Create a progress bar window
    % Loop through each tetrahedron
    for i = 1:size(tetrahedra_with_remaining_edges, 1)
        v = tetrahedra_with_remaining_edges(i, :);
        % Update waitbar
        waitbar(i / size(tetrahedra_with_remaining_edges, 1), holder, sprintf('Progress RBC: %d%%',round((i / size(tetrahedra_with_remaining_edges, 1)) * 100)));
        % Compute centroid (center point) of the tetrahedron
        center_point_t = mean(Points(v, :), 1);
        distance_to_origin = sqrt(sum((center_point_t-center_point).^2));
        show_RBC = abs((distance_to_origin/concentration_factor)*randn()) <= max([a, b, c]);
        if show_RBC
            % Define constants for RBCs
            d = 8; % in meters
            br = 1; % in meters
            h = 2.12; % in meters
            % Calculate P, Q, and R for RBCs
            P = -(d^2/2) + (h^2/2) * ((d^2/br^2) - 1) - h^2/2 * ((d^2/br^2) - 1) * sqrt(1 - (br^2/h^2));
            Q = P * (d^2/br^2) + (br^2/4) * (d^4/br^4 - 1);
            R = -P * (d^2/4) - d^4/16;
            [x_rbc,y_rbc,z_rbc] = meshgrid(-10+center_point_t(1):0.5:10+ ...
                center_point_t(1),-10+center_point_t(2):0.5:10+center_point_t(2),...
                -10+center_point_t(3):0.5:10+center_point_t(3));
            cell_rotations = rand(1, 3) *pi;
            % Compute transformed coordinates for RBCs
            x_rbc_rot = x_rbc-center_point_t(1);
            y_rbc_rot = y_rbc-center_point_t(2);
            z_rbc_rot = z_rbc-center_point_t(3);
            % Apply rotations for RBCs
            x_rbc_temp = x_rbc_rot;
            x_rbc_rot = x_rbc_temp * cos(cell_rotations(1)) - z_rbc_rot * sin(cell_rotations(1));
            z_rbc_rot = x_rbc_temp * sin(cell_rotations(1)) + z_rbc_rot * cos(cell_rotations(1));
            y_rbc_temp = y_rbc_rot;
            y_rbc_rot = y_rbc_temp * cos(cell_rotations(2)) + z_rbc_rot * sin(cell_rotations(2));
            z_rbc_rot = -y_rbc_temp * sin(cell_rotations(2)) + z_rbc_rot * cos(cell_rotations(2));
            x_rbc_temp = x_rbc_rot;
            x_rbc_rot = x_rbc_temp * cos(cell_rotations(3)) - y_rbc_rot * sin(cell_rotations(3));
            y_rbc_rot = x_rbc_temp * sin(cell_rotations(3)) + y_rbc_rot * cos(cell_rotations(3));
            eq = ((x_rbc_rot).^2 + (y_rbc_rot).^2 + ...
                (z_rbc_rot).^2).^2 ...
                + P * ((x_rbc_rot).^2 + (y_rbc_rot).^2) ...
                + Q * (z_rbc_rot).^2 + R;
            % Initialize v_rbc to store presence
            v_rbc = zeros(size(x_rbc));
            % Mark points inside RBC
            v_rbc(eq <= 0) = 1;
            % Extract coordinates where v_rbc is 1
            rbc_x = x_rbc(v_rbc == 1);
            rbc_y = y_rbc(v_rbc == 1);
            rbc_z = z_rbc(v_rbc == 1);
            % Combine coordinates into a matrix
            rbc_coordinates = [rbc_x(:), rbc_y(:), rbc_z(:)];
            DTr = delaunayTriangulation(rbc_x, rbc_y, rbc_z);
            [Cr,vr] = convexHull(DTr);
            trisurf(Cr, DTr.Points(:,1), DTr.Points(:,2), DTr.Points(:,3), 'FaceColor', 'red', 'EdgeColor', 'none')
            % Generate material indices for each point
            material_index_rbc = ones(size(rbc_coordinates, 1), 1) * 1; % Material index for RBC
            % Append RBC coordinates to all_coordinates
            all_coordinates = [all_coordinates; rbc_coordinates, material_index_rbc];
        end
    end

    close(holder);

    % Randomly exclude some points
    exclude_percentage_platelet = rand()*0.5;
    num_points_to_exclude = round(exclude_percentage_platelet * size(Points, 1));
    indices_to_exclude = randperm(size(Points, 1), num_points_to_exclude);
    Points(indices_to_exclude, :) = [];



    platelet_radius = 1.5;
    for i = 1:size(Points,1)

        % Generate surface coordinates for the sphere
        [sx, sy, sz] = sphere;
        sx = sx * platelet_radius + Points(i, 1);
        sy = sy * platelet_radius + Points(i, 2);
        sz = sz * platelet_radius + Points(i, 3);
        surf(sx,sy,sz,'FaceColor', 'green', 'EdgeColor', 'none');
        pause(1e-12);
        % Append surface sphere coordinates to all_coordinates
        sphere_surface_coordinates = [sx(:), sy(:), sz(:), ones(size(sx(:)))*3]; % Material index for spheres
        all_coordinates = [all_coordinates; sphere_surface_coordinates]; %#ok<*AGROW>

        % Generate inside coordinates for the sphere
        num_inside_points = 500; % Adjust as needed
        sphere_center = Points(i, :);
        sphere_radius = platelet_radius;

        % Generate points uniformly distributed inside the sphere
        rand_points = rand(num_inside_points, 3) - 0.5; % Points in [-0.5, 0.5]
        rand_points = rand_points ./ sqrt(sum(rand_points.^2, 2)); % Normalize to unit vectors
        inside_sphere_coordinates = sphere_center + sphere_radius * rand_points;

        % Append inside sphere coordinates to all_coordinates
        sphere_inside_coordinates = [inside_sphere_coordinates, ones(size(inside_sphere_coordinates, 1), 1)*3]; % Material index for spheres
        all_coordinates = [all_coordinates; sphere_inside_coordinates];
    end

    % Sort the coordinates matrix based on first x, then y, and finally z
    all_coordinates_sorted = sortrows(single(all_coordinates), [4,1,2,3]);
    % % Save sorted coordinates to a text file
    writematrix(all_coordinates_sorted, sprintf('all_coordinates_%d.txt',index), 'Delimiter', 'tab');

    coordinates = all_coordinates_sorted(:,1:3);
    materials = all_coordinates_sorted(:,4);

    % safe distance from the boundaries in um 
    distance_to_boundary = 50; 

    % Define the size of the mesh space
    mesh_size_x = round(max(coordinates(:,1))-min(coordinates(:,1)))+2*distance_to_boundary;  %-+50um safe distance from boundary
    mesh_size_y = round(max(coordinates(:,2))-min(coordinates(:,2)))+2*distance_to_boundary;  %-+50um safe distance from boundary
    mesh_size_z = round(max(coordinates(:,3))-min(coordinates(:,3)))+2*distance_to_boundary;  %-+50um safe distance from boundary

    % Create an empty mesh space
    mesh_space = zeros(mesh_size_x, mesh_size_y, mesh_size_z, 'single');

    % Translate coordinates to positive quadrant
    coordinates = coordinates - min(coordinates) + 1;

    % Scale coordinates to fit within the mesh size
    scaling_factor_x = 1;
    scaling_factor_y = 1;
    scaling_factor_z = 1;
    coordinates_scaled(:,1) = coordinates(:,1) * scaling_factor_x + distance_to_boundary;
    coordinates_scaled(:,2) = coordinates(:,2) * scaling_factor_y + distance_to_boundary;
    coordinates_scaled(:,3) = coordinates(:,3) * scaling_factor_z + distance_to_boundary;

    % Ensure all coordinates are integers
    coordinates_scaled = round(coordinates_scaled);

    % Assign material indices to the mesh space
    for i = 1:size(coordinates_scaled, 1)
        x = coordinates_scaled(i, 1);
        y = coordinates_scaled(i, 2);
        z = coordinates_scaled(i, 3);

        % Check if the coordinates are within the mesh space
        if x >= 1 && x <= mesh_size_x && y >= 1 && y <= mesh_size_y && z >= 1 && z <= mesh_size_z
            % Assign the material index to the corresponding voxel
            mesh_space(x, y, z) = materials(i);
        end
    end

    % Save the mesh space to a file
    save(sprintf('mesh_space_%d.mat',index), 'mesh_space', '-v7.3');
    savefig(sprintf('clot_%d.fig',index));
end