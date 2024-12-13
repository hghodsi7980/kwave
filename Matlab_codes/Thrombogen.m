% ==========================================================
% Fibrin Network Simulation Script
% ==========================================================
function Thrombogen(counter)
seed = 10;
    fprintf("progress: %f\n",100*counter/18330);
    go_forward = false;
    seed = seed+counter;
    while ~ go_forward
        % Set the random seed for reproducibility
        seed = seed+1;
        rng(seed);

        % ==========================================================
        % Section 0: Initialize Parameters
        % ==========================================================
        numSpheres = max(randi(15),5); % Number of spherical inclusions
        max_diameter_factor = 0.5; % Maximum sphere diameter as a fraction of the largest dimension
        fibrin_concentration = 0.4+0.8*rand() ; % Generate random fibrin concentration between 0.2 and 2 g/L
        clot_volume = 1e6+2e7*rand(); % Generate a random clot volume between 8e3 and 3.7e7 um3
        Window_size = 150; % in um3
        platelet_ratio = 0.3; % Ratio of platelets to be plotted
        density_threshold = 0.001; % Minimum density threshold for nodes between inclusions (nodes/um^3)
        scaling_factor_RBC = 1;
        rbc_diameter = 8 / scaling_factor_RBC;  % Approximate RBC size
        % Set parameters based on fibrin concentration using interpolation (from literature data)
        C1 = 0.4; L1 = 4.87; Z1 = 3.19; B1_1 = 0.94; B_hat1 = 1.89; v1 = 0.339;
        C2 = 1.6; L2 = 2.99; Z2 = 3.33; B1_2 = 0.91; B_hat2 = 1.51; v2 = 0.341;
        % Define a prefered directional vector
        preferred_direction = [1, 1, 1];
        preferred_direction = preferred_direction / norm(preferred_direction); % Normalize
        direction_weight = 0; % Weight factor for directional alignment (adjust this based on your preference)

        % ==========================================================
        % Section 1: Generate Inclusions
        % ==========================================================
        % Generate spherical inclusions and add points around them
        max_dimension = clot_volume^(1/3);
        max_diameter = max_dimension * max_diameter_factor;
        sphere_centers = [];
        sphere_radii = [];
        Points = [];

        % Generate spheres iteratively, ensuring no overlap
        for i = 1:numSpheres
            while true
                % Calculate the average center of existing spheres, if any
                if isempty(sphere_centers)
                    avg_center = [0, 0, 0]; % Start at the origin if no spheres exist
                else
                    avg_center = mean(sphere_centers, 1);
                end

                % Generate random offset for sphere placement, relative to the average center
                distance_factor = rand(); % A random factor to control distance from the average center
                center_offset = (2 * rand(1, 3) - 1) * (1 - distance_factor) * (Window_size / 2);
                center = avg_center + center_offset;

                % Larger spheres closer to the average center, smaller further out
                if distance_factor < 0.9
                    radius = (0.5 + rand() * 0.5) * max_diameter; % Larger spheres (50%-100% of max diameter)
                else
                    radius = (0.1 + rand() * 0.4) * max_diameter; % Smaller spheres (10%-50% of max diameter)
                end

                % Check if this sphere intersects any previously added spheres
                intersects = false;
                for j = 1:length(sphere_radii)
                    distance = norm(center - sphere_centers(j, :));
                    if distance < (radius + sphere_radii(j))
                        intersects = true;
                        break;
                    end
                end

                % If no intersection, add the sphere to the list
                if ~intersects
                    sphere_centers = [sphere_centers; center];
                    sphere_radii = [sphere_radii; radius];
                    break;
                end
            end
        end
        % ==========================================================
        % Section 2: Calculate Effective Clot Volume
        % ==========================================================
        % Calculate total volume of spherical inclusions
        clot_volume_um3 = calculate_clot_volume(sphere_centers, sphere_radii);
        %disp(['Real clot volume :',num2str(clot_volume_um3)]);
        total_inclusion_volume_um3 = sum((4 / 3) * pi * (sphere_radii.^3));
        % Calculate effective clot volume (excluding inclusions)
        effective_clot_volume_um3 = clot_volume_um3 - total_inclusion_volume_um3;

        % Convert fibrin concentration from g/L to g/um^3 and estimate fibrin volume
        fibrin_concentration_um3 = fibrin_concentration * 1e-15; % Convert to g/um^3
        fibrin_density = 1.395 * 1e-12; % Assume fibrin density in g/um^3
        fibrin_volume_um3 = (fibrin_concentration_um3 / fibrin_density) * effective_clot_volume_um3;

        % ==========================================================
        % Section 3: Estimate Bond Length and Number of Points
        % ==========================================================
        % Calculate total bond length required (assuming cylindrical bonds)
        fibrin_radius = 0.135 / 2; % um

        % Total bond length required
        total_bond_length = fibrin_volume_um3 / (pi * fibrin_radius^2);

        % Calculate target mean bond length and estimate the number of required nodes
        target_mean_length = interp1([C1, C2], [L1, L2], fibrin_concentration, 'linear', 'extrap');
        numPoints = ceil(total_bond_length / target_mean_length);
        if numPoints < 1000
            go_forward = false;
        end
        % ==========================================================
        % Section 4: Generate Cluster Points Around Inclusions
        % ==========================================================
        % Generate points in clusters around spherical inclusions
        for i = 1:numSpheres
            cluster_center = sphere_centers(i, :);
            radius = sphere_radii(i);
            num_cluster_points = round(numPoints / numSpheres);
            count = 0;
            cluster_points = [];
            bandwidth_std = eye(3) * (clot_volume_um3 / numSpheres) / 250;  % Ensure connectivity
            while count < num_cluster_points
                % Generate points using Gaussian distribution
                new_points = mvnrnd(cluster_center, bandwidth_std, num_cluster_points);

                % Calculate distances of points from the current sphere center
                distances = sqrt(sum((new_points - cluster_center).^2, 2));

                % Keep only points that are outside the current sphere radius
                valid_points = new_points(distances > radius, :);

                % Ensure valid points are not inside any other spheres
                for q = 1:length(sphere_radii)
                    if q ~= i
                        other_sphere_center = sphere_centers(q, :);
                        other_sphere_radius = sphere_radii(q);
                        distances_to_other_sphere = sqrt(sum((valid_points - other_sphere_center).^2, 2));
                        valid_points = valid_points(distances_to_other_sphere > other_sphere_radius, :);
                    end
                end

                % Append valid points until reaching the desired number of points
                points_needed = min(num_cluster_points - count, size(valid_points, 1));
                cluster_points = [cluster_points; valid_points(1:points_needed, :)];
                count = count + points_needed;
            end

            % Append to overall points list
            Points = [Points; cluster_points];
        end

        density = plot_density_based_points(Points, 10);
        if ~isempty(density)
            std_density = std(density(:,2)/max(density(:,2)));
        else
            std_density = 0;
        end

        if std_density > 0.25
            %disp('Not good enough!')
            %disp(std_density)
            go_forward = false;
        else
            go_forward = true;
        end
    end

    % ==========================================================
    % Section 5: Generate Bonds Based on Valency Distribution
    % ==========================================================
    % Generate bonds between points to achieve the target valency
    target_valency = interp1([C1, C2], [Z1, Z2], fibrin_concentration, 'linear', 'extrap');
    [bonds, bond_lengths] = generate_bonds_with_valency(Points, target_valency, preferred_direction,direction_weight);

    % Calculate the number of connections for each point
    numPoints = size(Points, 1); % Update to reflect the current number of points
    connection_counts = zeros(numPoints, 1); % Initialize connection counts
    for i = 1:size(bonds, 1)
        connection_counts(bonds(i, 1)) = connection_counts(bonds(i, 1)) + 1;
        connection_counts(bonds(i, 2)) = connection_counts(bonds(i, 2)) + 1;
    end

    % ==========================================================
    % Section 6: Squeeze Network Based on Sphere Inclusions
    % ==========================================================
    for i = 1:size(Points, 1)
        for j = 1:length(sphere_radii)
            dist_to_center = norm(Points(i, :) - sphere_centers(j, :));
            if dist_to_center < sphere_radii(j) * 1.2  % 20% larger than radius to create an influence zone
                % Move the point outwards away from the sphere center
                direction = (Points(i, :) - sphere_centers(j, :)) / dist_to_center;
                displacement = (sphere_radii(j) * 1.2 - dist_to_center) * direction;
                Points(i, :) = Points(i, :) + displacement;
            end
        end
    end

    % ==========================================================
    % Section 7: Strain Energy Relaxation with Balanced Penalties
    % ==========================================================
    % Parameters for strain energy relaxation and optimization
    max_iterations = 1000;
    energy_threshold = 0.07;
    k_fibrin = 1;
    initial_step_size = 0.001;
    step_size = initial_step_size;
    max_step_size = 0.05;
    min_step_size = 1e-5;

    % New parameters for convergence stability
    kl_stable_threshold = 1e-2;  % Threshold for considering J-S Divergence as stable
    kl_stable_max_count = 5;    % Number of successive iterations with small changes to consider convergence stable

    % Initialize previous J-S Divergence for comparison in the iteration loop
    prev_Jensen_Shannon_Divergence = Inf; % Start with a high initial value

    % Pre-calculate reference bond length distribution for comparison
    nu = interp1([C1, C2], [v1, v2], fibrin_concentration, 'linear', 'extrap');
    L = target_mean_length; % Mean length
    % Calculate s^2 from nu
    s_squared = nu * L^2;
    % Calculate zeta^2
    zeta_squared = log((s_squared / L^2) + 1);
    % Calculate lambda
    lambda = log(L) - (zeta_squared / 2);
    % Define the range of x values for the distribution
    bond_lengths = compute_bond_lengths(bonds, Points);
    x_values = linspace(1e-10, max(bond_lengths), length(bond_lengths)); % Avoid zero for stability
    % Calculate the log-normal PDF based on the paper formula
    log_normal_pdf = (1 ./ (x_values * sqrt(2 * pi * zeta_squared))) .* ...
        exp(-((lambda - log(x_values)).^2) / (2 * zeta_squared));

    % ==========================================================
    % Modified Relaxation Process to Include Directional Force
    % ==========================================================
    for iteration = 1:max_iterations
        % Compute current bond lengths
        bond_lengths = compute_bond_lengths(bonds, Points);

        % Fit Log-normal distributions and compute PDFs
        pd1 = fitdist(bond_lengths, 'Lognormal');
        dx = mean(diff(x_values));
        pdf1 = pdf(pd1, x_values);

        % Normalize PDFs
        pdf1_normalized = pdf1 / sum(pdf1 * dx);
        log_normal_pdf_normalized = log_normal_pdf / sum(log_normal_pdf * dx);

        % Calculate average distribution
        average_pdf = 0.5 * (pdf1_normalized + log_normal_pdf_normalized);

        % Add epsilon to avoid log(0)
        epsilon = 1e-10;
        pdf1_normalized = pdf1_normalized + epsilon;
        log_normal_pdf_normalized = log_normal_pdf_normalized + epsilon;
        average_pdf = average_pdf + epsilon;

        % Compute J-S Divergence
        KL_P1_avg = sum(pdf1_normalized .* log(pdf1_normalized ./ average_pdf)) * dx;
        KL_P2_avg = sum(log_normal_pdf_normalized .* log(log_normal_pdf_normalized ./ average_pdf)) * dx;
        Jensen_Shannon_Divergence = 0.5 * (KL_P1_avg + KL_P2_avg);

        if isnan(Jensen_Shannon_Divergence)
            %disp('NaN detected in Jensen-Shannon Divergence. Stopping iteration.');
            break;
        end

        %disp(['Jensen-Shannon Divergence: ', num2str(Jensen_Shannon_Divergence)]);

        % Compute forces considering the current bond length distribution and move points
        forces = compute_forces_log_normal(bonds, bond_lengths, lambda, zeta_squared, k_fibrin, Points);

        % Add directional alignment force
        directional_forces = zeros(size(forces));
        for i = 1:size(bonds, 1)
            p1 = bonds(i, 1);
            p2 = bonds(i, 2);
            bond_vector = Points(p2, :) - Points(p1, :);
            bond_vector_normalized = bond_vector / norm(bond_vector);

            % Calculate directional force as a projection onto the preferred direction
            alignment_component = dot(bond_vector_normalized, preferred_direction);
            alignment_force = direction_weight * alignment_component * preferred_direction;

            % Apply directional force
            directional_forces(p1, :) = directional_forces(p1, :) + alignment_force;
            directional_forces(p2, :) = directional_forces(p2, :) - alignment_force;
        end

        % Total force = original force + directional force
        total_forces = forces + directional_forces;

        % Move points based on computed total forces
        Points = move_points(Points, total_forces, step_size);

        % Adjust step size and k_fibrin
        if iteration > 1 && Jensen_Shannon_Divergence >= prev_Jensen_Shannon_Divergence
            step_size = max(step_size * 0.8, min_step_size);
            k_fibrin = min(k_fibrin * 2, 20); % Increase k_fibrin significantly for more exploration
            %disp(['Step size reduced to ', num2str(step_size), ', k_fibrin increased to ', num2str(k_fibrin), ' at iteration ', num2str(iteration)]);
        else
            step_size = min(step_size * 1.1, max_step_size);
            k_fibrin = max(k_fibrin * 0.9, 1); % Gradually reduce k_fibrin if improving
        end

        % Check for convergence and update kl_stable_count
        if abs(Jensen_Shannon_Divergence - prev_Jensen_Shannon_Divergence) < kl_stable_threshold
            kl_stable_count = kl_stable_count + 1;
        else
            kl_stable_count = 0;
        end

        % Apply random perturbation if the system is stable for too many iterations
        if kl_stable_count > kl_stable_max_count
            % Apply a small random perturbation to all points
            perturbation_strength = 0.05 * target_mean_length; % Relative to target mean length
            Points = Points + perturbation_strength * (rand(size(Points)) - 0.5);
            %disp(['Applied random perturbations to points at iteration ', num2str(iteration)]);
            kl_stable_count = 0; % Reset the stable count
        end


        % Update Jensen-Shannon Divergence for the next iteration
        prev_Jensen_Shannon_Divergence = Jensen_Shannon_Divergence;

        % Check if the system has converged
        if Jensen_Shannon_Divergence < energy_threshold
            %disp(['Converged at iteration ', num2str(iteration)]);
            break;
        end
    end


    % ==========================================================
    % Section 8: Post-Optimization Cleanup for Long Bonds
    % ==========================================================
    max_allowed_length = 5 * target_mean_length; % Maximum allowed bond length (adjust as needed)

    % Identify and remove nodes with many super-long bonds
    long_bond_threshold = max_allowed_length;
    to_remove = false(numPoints, 1); % Logical array to mark nodes for removal

    % Loop over each node to identify problematic ones
    for i = 1:numPoints
        % Find bonds connected to the current node
        connected_bonds = bonds(:, 1) == i | bonds(:, 2) == i;
        connected_lengths = bond_lengths(connected_bonds);

        % Check if there are multiple super-long bonds connected to this node
        if sum(connected_lengths > long_bond_threshold) > 0
            % Mark the node for removal
            to_remove(i) = true;
        end
    end

    % Remove marked nodes from Points and update bonds
    Points(to_remove, :) = [];
    remaining_indices = find(~to_remove);

    % Update bonds to reflect the removal of nodes
    new_bonds = [];
    for i = 1:size(bonds, 1)
        % Get the bond nodes and check if both are still valid
        node1 = bonds(i, 1);
        node2 = bonds(i, 2);
        if ~to_remove(node1) && ~to_remove(node2)
            % Map old indices to new indices after removal
            new_node1 = find(remaining_indices == node1);
            new_node2 = find(remaining_indices == node2);
            new_bonds = [new_bonds; new_node1, new_node2];
        end
    end

    % Update bonds and recalculate bond lengths
    bonds = new_bonds;
    bond_lengths = compute_bond_lengths(bonds, Points);

    % ==========================================================
    % Section: Remove Disconnected Components
    % ==========================================================
    % Create an adjacency list representation of the graph
    numPoints = size(Points, 1);
    adjacency_list = cell(numPoints, 1);
    for i = 1:size(bonds, 1)
        adjacency_list{bonds(i, 1)} = [adjacency_list{bonds(i, 1)}, bonds(i, 2)];
        adjacency_list{bonds(i, 2)} = [adjacency_list{bonds(i, 2)}, bonds(i, 1)];
    end

    % Find connected components using BFS
    visited = false(numPoints, 1);
    components = {};
    for i = 1:numPoints
        if ~visited(i)
            % Start a new BFS from the unvisited node
            queue = i;
            component = [];
            while ~isempty(queue)
                current = queue(1);
                queue(1) = [];
                if ~visited(current)
                    visited(current) = true;
                    component = [component, current];
                    queue = [queue, adjacency_list{current}]; % Add all neighbors to the queue
                end
            end
            components{end + 1} = component; % Store the connected component
        end
    end

    % Find the largest connected component
    component_sizes = cellfun(@length, components);
    [~, largest_component_idx] = max(component_sizes);
    largest_component = components{largest_component_idx};

    % Keep only the points and bonds that are part of the largest connected component
    Points = Points(largest_component, :);

    % Update the bonds to only include bonds within the largest connected component
    % Create a mapping from old indices to new indices
    new_index_map = zeros(numPoints, 1);
    new_index_map(largest_component) = 1:length(largest_component);

    % Update the bonds array
    bonds = bonds(ismember(bonds(:, 1), largest_component) & ismember(bonds(:, 2), largest_component), :);
    bonds = [new_index_map(bonds(:, 1)), new_index_map(bonds(:, 2))];

    % ==========================================================
    % Section 10: move the Bonds Inside Spheres
    % ==========================================================
    % move bonds that lie inside the spherical inclusions
    bonds = move_bonds_outside_spheres(bonds, Points, sphere_centers, sphere_radii);

    % Calculate the number of connections for each point
    numPoints = size(Points, 1); % Update to reflect the current number of points
    connection_counts = zeros(numPoints, 1); % Initialize connection countsg
    for i = 1:size(bonds, 1)
        connection_counts(bonds(i, 1)) = connection_counts(bonds(i, 1)) + 1;
        connection_counts(bonds(i, 2)) = connection_counts(bonds(i, 2)) + 1;
    end

    % ==========================================================
    % Section 11: Plot Final Network with Bezier Curves, Platelets, and RBCs in Sphere Inclusions
    % ==========================================================
    for i = 1:size(bonds, 1)
        P0 = Points(bonds(i, 1), :);  % Start point of the bond
        P2 = Points(bonds(i, 2), :);  % End point of the bond
        midpoint = (P0 + P2) / 2;  % Midpoint for control point
        random_curvature = (rand(1, 3) - 0.5)*2;  % Random curvature
        P1 = midpoint + random_curvature;  % Control point for Bezier curve
        num_bond_points = 20;
        t = linspace(0, 1, num_bond_points);  % Parameter for Bezier curve
        bezier_curve = (1 - t').^2 * P0 + 2 * (1 - t') .* t' * P1 + t'.^2 * P2;  % Bezier curve formula
    end
    % Plot platelets at nodes with 4 connections (cross-links)
    platelet_indices = find(connection_counts == 4);
    num_platelets_to_plot = round(platelet_ratio * length(platelet_indices));
    selected_indices = randsample(platelet_indices, num_platelets_to_plot);

    min_distance = rbc_diameter;  % Distance between RBCs (just beside, no intersection)
    rbc_points_all = [];  % to save the RBC points
    for i = 1:length(sphere_radii)
        radius_sphere = sphere_radii(i);
        sphere_center = sphere_centers(i, :);
        num_RBCs = (radius_sphere/(rbc_diameter/2))^3/3 ;  % Number of RBCs
        % Step 3: Randomly generate positions inside the sphere and allow adjacent placement
        rbc_positions = [];
        % Place the first RBC randomly inside the sphere
        theta = 2 * pi * rand();
        phi = acos(2 * rand() - 1);
        r = radius_sphere * rand()^(1/3);  % Cubic root for uniform distribution inside sphere
        rbc_center = sphere_center + [r * sin(phi) * cos(theta), r * sin(phi) * sin(theta), r * cos(phi)];
        rbc_positions = [rbc_positions; rbc_center];

        % Place remaining RBCs adjacent to existing ones
        for j = 2:floor(num_RBCs)-1
            while true
                % Choose an existing RBC to connect with
                parent_rbc_idx = randi(size(rbc_positions, 1));
                parent_rbc_center = rbc_positions(parent_rbc_idx, :);

                % Generate a new point near the selected RBC
                offset_theta = 2 * pi * rand();
                offset_phi = acos(2 * rand() - 1);
                offset_distance = rbc_diameter + rand() * (rbc_diameter / 2);  % Small offset
                new_rbc_center = parent_rbc_center + ...
                    [offset_distance * sin(offset_phi) * cos(offset_theta), ...
                    offset_distance * sin(offset_phi) * sin(offset_theta), ...
                    offset_distance * cos(offset_phi)];

                % Ensure the new RBC is still within the sphere
                if norm(new_rbc_center - sphere_center) <= radius_sphere
                    % Ensure the new RBC doesn't intersect with any existing RBCs
                    distances = sqrt(sum((rbc_positions - new_rbc_center).^2, 2));
                    if all(distances > min_distance)
                        rbc_positions = [rbc_positions; new_rbc_center];  % Add new RBC position
                        break;
                    end
                end
            end
        end

        % Plot RBCs with biconcave shape in the current sphere inclusion
        for j = 1:size(rbc_positions, 1)
            center_point_t = rbc_positions(j, :);
            cell_rotations = rand(1, 3) * pi;

            % RBC code as provided (adapt to each position and rotation)
            d = 8 / scaling_factor_RBC;
            br = 1 / scaling_factor_RBC;
            h = 2.12 / scaling_factor_RBC;
            P = -(d^2 / 2) + (h^2 / 2) * ((d^2 / br^2) - 1) - h^2 / 2 * ((d^2 / br^2) - 1) * sqrt(1 - (br^2 / h^2));
            Q = P * (d^2 / br^2) + (br^2 / 4) * (d^4 / br^4 - 1);
            R = -P * (d^2 / 4) - d^4 / 16;

            [x_rbc, y_rbc, z_rbc] = meshgrid(-10 + center_point_t(1):0.5:10 + center_point_t(1), ...
                -10 + center_point_t(2):0.5:10 + center_point_t(2), ...
                -10 + center_point_t(3):0.5:10 + center_point_t(3));
            % Apply rotations
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

            % Equation for RBC shape
            eq = ((x_rbc_rot).^2 + (y_rbc_rot).^2 + (z_rbc_rot).^2).^2 + P * ((x_rbc_rot).^2 + (y_rbc_rot).^2) + Q * (z_rbc_rot).^2 + R;
            v_rbc = zeros(size(x_rbc));
            v_rbc(eq <= 0) = 1;
            rbc_x = x_rbc(v_rbc == 1);
            rbc_y = y_rbc(v_rbc == 1);
            rbc_z = z_rbc(v_rbc == 1);

            % Convex Hull for RBC shape
            rbc_points = [rbc_x(:), rbc_y(:), rbc_z(:)];
            if size(rbc_points, 1) >= 4  % Ensure there are enough points for convex hull
                K = convhull(rbc_points);  % Get convex hull
                rbc_points_all = [rbc_points_all ; rbc_points];
            end
        end
    end

    % ==========================================================
    % Define the cropping window and matrix size to save the clot block
    % ==========================================================
    Window_size_crop = 120; % Total range in each direction is 120um (-60 to 60)
    Resolution = 240; % Number of voxels per axis
    VoxelSize = Window_size_crop / Resolution; % Size of each voxel in um
    Window_center = mean(Points,1);

    % Define the bounds of the cropping window
    crop_min = Window_center-Window_size_crop/2; % Minimum bound in um
    crop_max = Window_center+Window_size_crop/2;  % Maximum bound in um

    % Initialize the 3D matrix
    ClotMatrix = zeros(Resolution, Resolution, Resolution);

    % Function to map coordinates into the matrix indices
    map_to_matrix = @(coords) round((coords - crop_min) / VoxelSize) + 1;

    % ==========================================================
    % Process fibrin bonds
    % ==========================================================
    for i = 1:size(bonds, 1)
        P1 = Points(bonds(i, 1), :);
        P2 = Points(bonds(i, 2), :);

        % Check if the bond is within the cropping window
        if all(P1 >= crop_min & P1 <= crop_max) && all(P2 >= crop_min & P2 <= crop_max)
            % Map bond points to matrix indices
            P1_idx = map_to_matrix(P1);
            P2_idx = map_to_matrix(P2);

            % Use Bresenham's algorithm to interpolate points along the bond
            bond_voxels = bresenham3d(P1_idx, P2_idx);
            for j = 1:size(bond_voxels, 1)
                if all(bond_voxels(j, :) > 0 & bond_voxels(j, :) <= Resolution)
                    ClotMatrix(bond_voxels(j, 1), bond_voxels(j, 2), bond_voxels(j, 3)) = 2; % Fibrin bond
                end
            end
        end
    end

    % ==========================================================
    % Process RBCs
    % ==========================================================
    for i = 1:size(rbc_points_all, 1)
        RBC = rbc_points_all(i, :);

        % Check if the RBC is within the cropping window
        if all(RBC >= crop_min & RBC <= crop_max)
            % Map RBC position to matrix indices
            RBC_idx = map_to_matrix(RBC);
            if all(RBC_idx > 0 & RBC_idx <= Resolution)
                ClotMatrix(RBC_idx(1), RBC_idx(2), RBC_idx(3)) = 1; % RBC
            end
        end
    end

% ==========================================================
% Process platelets
% ==========================================================
platelet_radius = 1.5; % Radius of the platelet in microns

for i = 1:length(platelet_indices)
    platelet_pos = Points(platelet_indices(i), :);
    
    % Check if the platelet is within the cropping window
    if all(platelet_pos >= crop_min & platelet_pos <= crop_max)
        % Map platelet center position to matrix indices
        platelet_center_idx = map_to_matrix(platelet_pos);

        if all(platelet_center_idx > 0 & platelet_center_idx <= Resolution)
            % Generate a sphere of platelet points around the center
            [x_sphere, y_sphere, z_sphere] = ndgrid(...
                (platelet_pos(1) - platelet_radius):VoxelSize:(platelet_pos(1) + platelet_radius), ...
                (platelet_pos(2) - platelet_radius):VoxelSize:(platelet_pos(2) + platelet_radius), ...
                (platelet_pos(3) - platelet_radius):VoxelSize:(platelet_pos(3) + platelet_radius));

            % Flatten the grids and compute distances from the platelet center
            sphere_coords = [x_sphere(:), y_sphere(:), z_sphere(:)];
            distances = sqrt(sum((sphere_coords - platelet_pos).^2, 2));

            % Select points within the sphere radius
            valid_sphere_points = sphere_coords(distances <= platelet_radius, :);

            % Map valid points to matrix indices and update the ClotMatrix
            for j = 1:size(valid_sphere_points, 1)
                sphere_point_idx = map_to_matrix(valid_sphere_points(j, :));
                if all(sphere_point_idx > 0 & sphere_point_idx <= Resolution)
                    ClotMatrix(sphere_point_idx(1), sphere_point_idx(2), sphere_point_idx(3)) = 3; % Platelet
                end
            end
        end
    end
end

    % ==========================================================
    % Save the cropped matrix
    % ==========================================================
    save(sprintf('CroppedClotMatrix_%d.mat',counter), 'ClotMatrix');
end
% ==========================================================
% Section 12: Helper Functions
% ==========================================================
% Note: Helper functions like generate_bonds_with_valency, compute_forces_log_normal,
% compute_bond_lengths, move_points, remove_bonds_inside_spheres, and compute_kl_divergence
% are defined below to perform specific operations.
% They are intended to modularize the workflow and make the main code easier to read and maintain.

% ==========================================================
% Function: generate_bonds_with_valency (modified for directional preference)
% ==========================================================
function [bonds, bond_lengths] = generate_bonds_with_valency(Points, target_valency, preferred_direction, direction_weight)
numPoints = size(Points, 1);
bonds = [];
connection_counts = zeros(numPoints, 1); % Track the number of connections per point
target_connections = round(geornd(1 /(target_valency-2), numPoints, 1)+target_valency); % Shifted geometric distribution

% Generate bonds to match the target valency distribution
for i = 1:numPoints
    % Calculate distances from the current point to all other points
    distances = sqrt(sum((Points - Points(i, :)).^2, 2));

    % Calculate the dot product with preferred direction to favor alignment
    vectors = Points - Points(i, :);
    direction_scores = vectors * preferred_direction'; % Projection onto preferred direction

    % Create a weighted distance metric
    weighted_distances = distances - direction_weight * direction_scores;

    [~, idx] = sort(weighted_distances);
    nearest_neighbors = idx(2:end); % Nearest neighbors, excluding itself

    % Connect to nearest neighbors until reaching the target valency
    for j = 1:length(nearest_neighbors)
        if connection_counts(i) < target_connections(i) && connection_counts(nearest_neighbors(j)) < target_connections(i)
            bonds = [bonds; i, nearest_neighbors(j)];
            connection_counts(i) = connection_counts(i) + 1;
            connection_counts(nearest_neighbors(j)) = connection_counts(nearest_neighbors(j)) + 1;
        end
        % Stop if the target valency is reached
        if connection_counts(i) >= target_connections(i)
            break;
        end
    end
end

% Calculate bond lengths for each generated bond
bond_lengths = compute_bond_lengths(bonds, Points);
end





% ==========================================================
% Function: compute_forces_log_normal
% ==========================================================
function forces = compute_forces_log_normal(bonds, bond_lengths, lambda, zeta_squared, k_fibrin, Points)
numPoints = size(Points, 1);
forces = zeros(numPoints, 3);

% Adjust the target length dynamically using a weighted average approach
target_length = exp(lambda + (zeta_squared / 2));
target_length_adjustment = 0.5 * mean(bond_lengths); % Dynamic adjustment to balance force application

for i = 1:size(bonds, 1)
    p1 = bonds(i, 1);
    p2 = bonds(i, 2);
    bond_vector = Points(p2, :) - Points(p1, :);
    bond_length = bond_lengths(i);

    % Dynamic target to ensure balance between compression and expansion
    adjusted_target_length = 0.5 * (target_length + target_length_adjustment);

    % Force calculation to drive towards target length
    force_magnitude = -k_fibrin * (bond_length - adjusted_target_length);
    force_vector = (force_magnitude / bond_length) * bond_vector;

    % Apply forces
    forces(p1, :) = forces(p1, :) - force_vector;
    forces(p2, :) = forces(p2, :) + force_vector;
end
end

% ==========================================================
% Function: compute_bond_lengths
% ==========================================================
function bond_lengths = compute_bond_lengths(bonds, Points)
% Compute the lengths of bonds between points
epsilon = 1e-6; % Small value to ensure numerical stability
bond_lengths = zeros(size(bonds, 1), 1);
for i = 1:size(bonds, 1)
    bond_lengths(i) = norm(Points(bonds(i, 1), :) - Points(bonds(i, 2), :)) + epsilon;
end
end

% ==========================================================
% Function: move_points
% ==========================================================
function Points = move_points(Points, forces, step_size)
% Update point positions based on forces and step size
Points = Points + step_size * forces;
end

% ==========================================================
% Function: move_bonds_outside_spheres
% ==========================================================
function bonds = move_bonds_outside_spheres(bonds, Points, sphere_centers, sphere_radii)
% Move bonds that pass through spherical inclusions to ensure they are outside
numBonds = size(bonds, 1);
for i = 1:numBonds
    midpoint = (Points(bonds(i, 1), :) + Points(bonds(i, 2), :)) / 2;
    for j = 1:length(sphere_radii)
        if norm(midpoint - sphere_centers(j, :)) < sphere_radii(j)
            % Move midpoint out of the sphere
            direction = (midpoint - sphere_centers(j, :)) / norm(midpoint - sphere_centers(j, :));
            new_midpoint = sphere_centers(j, :) + direction * sphere_radii(j) * 1.1; % Move slightly outside the sphere radius

            % Update the position of the midpoint
            Points(bonds(i, 1), :) = (Points(bonds(i, 1), :) + new_midpoint) / 2;
            Points(bonds(i, 2), :) = (Points(bonds(i, 2), :) + new_midpoint) / 2;
        end
    end
end
end

% ==========================================================
% Function: calculate_local_density
% ==========================================================
function [point_density] = plot_density_based_points(Points, radius)
% Plot the points with colors based on their local density
% Inputs:
%   Points - an Nx3 array of point coordinates
%   radius - radius for determining local density
%   volume_size - size of the overall cubic volume for normalization purposes
% Outputs:
%   point_density - a 2xN array representing each point's coordinates and its local density

numPoints = size(Points, 1);
densities = zeros(numPoints, 1);

% Calculate the local density for each point
for i = 1:numPoints
    % Calculate distances to all other points
    distances = sqrt(sum((Points - Points(i, :)).^2, 2));

    % Count how many points are within the specified radius, excluding the point itself
    point_count = sum(distances < radius) - 1;

    % Calculate the density as the number of points within the radius divided by the total number of points
    densities(i) = point_count / numPoints;
end

% Store density in the output matrix (point index and density value)
point_density = [transpose(1:numPoints), densities];

% Normalize densities for coloring purposes
densities_normalized = (densities - min(densities)) / (max(densities) - min(densities));

% Create a colormap for density-based coloring
cmap = jet(256); % Use the 'jet' colormap with 256 colors

end

% ==========================================================
% Function: calculate_clot volume
% ==========================================================
function total_volume = calculate_clot_volume(sphere_centers, sphere_radii)
total_points = 1000;
% Check if sphere_centers is N by 3 and sphere_radii is N by 1
if size(sphere_centers, 2) ~= 3 || size(sphere_centers, 1) ~= size(sphere_radii, 1)
    error('The number of sphere centers and radii must match, and centers must be a N by 3 matrix.');
end

num_spheres = length(sphere_radii);

% Calculate the surface area of each sphere
surface_areas = 4 * pi * (sphere_radii.^2);

% Distribute the points proportionally based on surface area
points_per_sphere = round((surface_areas / sum(surface_areas)) * total_points);

all_points = [];

% Generate points on the surface of each sphere based on the allocated number
for i = 1:num_spheres
    radius = sphere_radii(i);
    center = sphere_centers(i, :);
    num_points = points_per_sphere(i);

    % Random points on a sphere surface using spherical coordinates
    theta = 2 * pi * rand(num_points, 1);     % Azimuthal angle
    phi = acos(2 * rand(num_points, 1) - 1);  % Polar angle

    x = radius * sin(phi) .* cos(theta) + center(1);
    y = radius * sin(phi) .* sin(theta) + center(2);
    z = radius * cos(phi) + center(3);

    % Append the points to the all_points matrix
    all_points = [all_points; [x, y, z]];
end

% Perform Delaunay Triangulation on all points
dt = delaunayTriangulation(all_points);

% Calculate the convex hull of the triangulated points
[K, volume] = convexHull(dt);

% Multiply the total volume by 1.1 as required
total_volume = volume * 1.1;
end

% ==========================================================
% Function: Bresenham3D Function
% ==========================================================

function line = bresenham3d(p1, p2)
% Bresenham-like algorithm for 3D line
% Inputs:
%   p1 - Start point [x, y, z]
%   p2 - End point [x, y, z]
% Output:
%   line - Nx3 array of voxel indices along the line

% Initialize variables
x1 = p1(1); y1 = p1(2); z1 = p1(3);
x2 = p2(1); y2 = p2(2); z2 = p2(3);
dx = abs(x2 - x1); dy = abs(y2 - y1); dz = abs(z2 - z1);
xs = sign(x2 - x1); ys = sign(y2 - y1); zs = sign(z2 - z1);

% Initialize line array
line = [];

% Determine dominant axis
if dx >= dy && dx >= dz
    % X-axis is dominant
    p1 = 2 * dy - dx;
    p2 = 2 * dz - dx;
    while x1 ~= x2
        line = [line; x1, y1, z1]; %#ok<AGROW>
        x1 = x1 + xs;
        if p1 >= 0
            y1 = y1 + ys;
            p1 = p1 - 2 * dx;
        end
        if p2 >= 0
            z1 = z1 + zs;
            p2 = p2 - 2 * dx;
        end
        p1 = p1 + 2 * dy;
        p2 = p2 + 2 * dz;
    end
elseif dy >= dx && dy >= dz
    % Y-axis is dominant
    p1 = 2 * dx - dy;
    p2 = 2 * dz - dy;
    while y1 ~= y2
        line = [line; x1, y1, z1]; %#ok<AGROW>
        y1 = y1 + ys;
        if p1 >= 0
            x1 = x1 + xs;
            p1 = p1 - 2 * dy;
        end
        if p2 >= 0
            z1 = z1 + zs;
            p2 = p2 - 2 * dy;
        end
        p1 = p1 + 2 * dx;
        p2 = p2 + 2 * dz;
    end
else
    % Z-axis is dominant
    p1 = 2 * dy - dz;
    p2 = 2 * dx - dz;
    while z1 ~= z2
        line = [line; x1, y1, z1]; %#ok<AGROW>
        z1 = z1 + zs;
        if p1 >= 0
            y1 = y1 + ys;
            p1 = p1 - 2 * dz;
        end
        if p2 >= 0
            x1 = x1 + xs;
            p2 = p2 - 2 * dz;
        end
        p1 = p1 + 2 * dy;
        p2 = p2 + 2 * dx;
    end
end

% Add the last point
line = [line; x2, y2, z2];
end



