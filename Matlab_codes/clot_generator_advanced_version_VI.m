clear; clc;close all
% -------------------
% Fixed random seed
% -------------------
seed = 50;
rng(seed);  % Set the random seed
% -------------------
% Parameters
% -------------------
% Generate numPoints with a normal distribution
numPoints = round(max(500, abs(normrnd(250, 1000)))); % Fibrin points
% Generate number of spheres (platelets) with a normal distribution
numSpheres = round(max(0, abs(normrnd(100, 40))));

max_iterations = 1000; % Max number of optimization iterations
energy_threshold = 1e-3; % Energy threshold for optimization stop
k_fibrin = 1; % Stiffness constant for fibrin fibers (in nN/µm)
step_size = 0.01; % Step size for moving points during optimization
max_diameter_factor = 0.5; % Maximum sphere diameter as a fraction of the largest dimension
platelet_diameter = 3; % Diameter of platelet spheres (3 um)
grid_size = 400; % Fixed size of 3D mesh



% ==========================================================
% Section 1: Initialization of Fibrin Network Using Random Clusters
% ==========================================================
% Generate numPoints with a normal distribution
numPoints = round(max(1000, abs(normrnd(5000, 2000)))); % Fibrin points
% Generate number of clusters
numClusters = randi([5, 10]); % Number of clusters
Points = [];

% Bonding box size 
a = round(max(20, abs(normrnd(75, 25))));
b = round(max(20, abs(normrnd(75, 25))));
c = round(max(20, abs(normrnd(75, 25))));

% Generate random clusters within the bounding box [a, b, c]
for i = 1:numClusters
    cluster_center = [a * rand(), b * rand(), c * rand()];
    % Generate points around each cluster center with multivariate normal distribution
    cluster_points = mvnrnd(cluster_center, eye(3) * 50, round(numPoints / numClusters));
    Points = [Points; cluster_points];
end

% Apply random displacements to further randomize the points
displacements = [randn(size(Points, 1), 1) * a / 10, randn(size(Points, 1), 1) * b / 10, randn(size(Points, 1), 1) * c / 10];
Points = Points + displacements;

% Calculate the clot volume using Delaunay triangulation and convex hull
DT = delaunayTriangulation(Points);
[C, clot_volume] = convexHull(DT);

disp(['Clot volume (from convex hull): ', num2str(clot_volume)]);

% Calculate the equilibrium bond length based on volume and number of points
d0_fibrin = (clot_volume / numPoints)^(1/3);
disp(['Calculated equilibrium bond length (d0_fibrin): ', num2str(d0_fibrin)]);


% ==========================================================
% Section 2: Strain Energy Relaxation of the Network
% ==========================================================
compute_energy = @(distances) sum(k_fibrin * (distances - d0_fibrin).^2);

% Generate initial bonds between points (Delaunay triangulation or nearest neighbors)
[bonds, bond_lengths] = generate_initial_bonds(Points);

for iteration = 1:max_iterations
    % Step 1: Compute current energy
    current_energy = compute_energy(bond_lengths);

    % Check if energy is NaN or below the threshold
    if isnan(current_energy)
        disp(['NaN encountered at iteration ', num2str(iteration), '. Optimization terminated.']);
        break;
    elseif current_energy < energy_threshold
        disp(['Optimization finished at iteration ', num2str(iteration), ' due to energy threshold.']);
        break;
    end

    % Step 2: Compute forces and move points
    forces = compute_forces(bonds, bond_lengths, d0_fibrin, k_fibrin, Points);
    Points = move_points(Points, forces, step_size);

    % Recompute bond lengths after moving points
    bond_lengths = compute_bond_lengths(bonds, Points);

    % Display current energy for debugging
    disp(['Iteration ', num2str(iteration), ', Energy: ', num2str(current_energy)]);

    % Plot the PDF of the bond lengths dynamically
    if mod(iteration, 10) == 0 || iteration == 1
        figure(1); % Create or access figure with ID 1
        histogram(bond_lengths, 'Normalization', 'pdf');
        xlabel('Bond Length');
        ylabel('Probability Density Function (PDF)');
        title(['PDF of Bond Lengths at Iteration ', num2str(iteration)]);

        % Set dynamic axis limits
        xlim([0, max(bond_lengths) + 5]); % Set x-axis limit based on the maximum bond length
        ylim([0, max(ylim) * 1.2]); % Set y-axis limit with some margin

        drawnow; % Update the plot dynamically
    end

end

% ==========================================================
% Section 3: Generate Spheres, Remove Extra Bonds/Points, and Squeeze Network
% ==========================================================
% Generate random spheres distributed around the network center
sphere_centers = [];
sphere_radii = [];

% Parameters for generating spheres around the network center
centroid = mean(Points);
std_dev_factor = 0.2; % Factor to control how clustered the spheres are around the network center
max_dimension = max([a, b, c]);
max_diameter = max_dimension * max_diameter_factor;

for i = 1:numSpheres
    while true
        % Generate random sphere center using normal distribution around the network centroid
        center = centroid + std_dev_factor * [std(Points(:,1)) * randn(), ...
            std(Points(:,2)) * randn(), ...
            std(Points(:,3)) * randn()];
        radius = rand() * max_diameter / 2;

        % Check if this sphere intersects any previously added spheres
        intersects = false;
        for j = 1:length(sphere_radii)
            distance = norm(center - sphere_centers(j, :));
            if distance < (radius + sphere_radii(j))
                intersects = true;
                break;
            end
        end

        % If it doesn't intersect, add it
        if ~intersects
            sphere_centers = [sphere_centers; center];
            sphere_radii = [sphere_radii; radius];
            break;
        end
    end
end

% Remove nodes inside spheres and squeeze network
inside_sphere = false(size(Points, 1), 1);  % Mark nodes inside spheres

for i = 1:size(Points, 1)
    for j = 1:length(sphere_radii)
        dist_to_sphere_center = norm(Points(i, :) - sphere_centers(j, :));
        if dist_to_sphere_center < sphere_radii(j)
            inside_sphere(i) = true;
            break;
        end
    end
end

% Remove nodes inside spheres
Points(inside_sphere, :) = [];

% Update bond indices after removing nodes
new_indices = find(~inside_sphere);  % Indices of the remaining nodes
mapping = zeros(size(inside_sphere));  % Create a mapping from old to new indices
mapping(new_indices) = 1:length(new_indices);  % Fill in the mapping

% Update bonds based on the new mapping
bonds = bonds(~ismember(bonds(:,1), find(inside_sphere)) & ~ismember(bonds(:,2), find(inside_sphere)), :);
bonds = [mapping(bonds(:, 1)), mapping(bonds(:, 2))];  % Re-map bond indices to the new Points array

% Squeeze the network based on the sphere influence
for i = 1:size(Points, 1)
    for j = 1:length(sphere_radii)
        % Calculate distance from point to sphere center
        dist_to_center = norm(Points(i, :) - sphere_centers(j, :));
        if dist_to_center < sphere_radii(j) * 2  % Influence zone
            direction = (Points(i, :) - sphere_centers(j, :)) / dist_to_center;
            influence_factor = (sphere_radii(j) * 2 - dist_to_center) / (sphere_radii(j) * 2);
            Points(i, :) = Points(i, :) + influence_factor * direction * 5;  % Move the point outward
        end
    end
end

bond_lengths = compute_bond_lengths(bonds, Points);

figure;
histogram(bond_lengths, 'Normalization', 'pdf');
xlabel('Bond Length');
ylabel('Probability Density Function (PDF)');
title('PDF of Bond Lengths after squeezing');

% Set dynamic axis limits
xlim([0, max(bond_lengths) + 5]); % Set x-axis limit based on the maximum bond length
ylim([0, max(ylim) * 1.2]); % Set y-axis limit with some margin


% ==========================================================
% Section 4: Generate Platelets at Nodes with 4 Connections
% ==========================================================
platelet_centers = [];

% Identify nodes with exactly 4 connections
connection_counts = zeros(size(Points, 1), 1);
for i = 1:size(bonds, 1)
    connection_counts(bonds(i, 1)) = connection_counts(bonds(i, 1)) + 1;
    connection_counts(bonds(i, 2)) = connection_counts(bonds(i, 2)) + 1;
end

% Find indices of nodes with exactly 4 connections
nodes_with_4_connections = find(connection_counts == 4);

% Place platelets at these nodes
for i = 1:length(nodes_with_4_connections)
    platelet_centers = [platelet_centers; Points(nodes_with_4_connections(i), :)];
end

% Display the number of platelets generated
disp(['Number of platelets generated: ', num2str(length(platelet_centers))]);


% ==========================================================
% Section 5: Generate Bezier Curve Points on Fibrin Bonds and Inside Platelets
% ==========================================================
num_bond_points = 50; % Number of points per fibrin bond (on Bezier curve)
num_platelet_points = 100; % Number of points on and inside each platelet

coordinate_matrix = [];

% Generate points on Bezier curves for each bond
for i = 1:size(bonds, 1)
    P0 = Points(bonds(i, 1), :);  % Start point of the bond
    P2 = Points(bonds(i, 2), :);  % End point of the bond
    midpoint = (P0 + P2) / 2;  % Midpoint for generating control point

    % Generate random curvature for Bezier curve
    random_curvature = (rand(1, 3) - 0.5) * 3;  % Small curvature
    P1 = midpoint + random_curvature;  % Control point

    % Generate points on the Bezier curve
    t = linspace(0, 1, num_bond_points);  % Parameter for Bezier curve
    bezier_points = (1 - t').^2 * P0 + 2 * (1 - t') .* t' * P1 + t'.^2 * P2;  % Bezier curve formula

    atom_type = ones(num_bond_points, 1);  % 1 for Fibrin points on the bond
    coordinate_matrix = [coordinate_matrix; [bezier_points, atom_type]];  % Append Bezier points to coordinate matrix
end

% Generate 50 points on and inside each platelet
for i = 1:size(platelet_centers, 1)
    center = platelet_centers(i, :);
    for j = 1:num_platelet_points
        % Random point inside the platelet sphere
        r = (rand()^(1/3)) * (platelet_diameter / 2); % Random distance within the sphere
        theta = 2 * pi * rand(); % Random azimuthal angle
        phi = acos(2 * rand() - 1); % Random polar angle
        point = center + [r * sin(phi) * cos(theta), r * sin(phi) * sin(theta), r * cos(phi)];
        atom_type = 2; % 2 for Platelet
        coordinate_matrix = [coordinate_matrix; [point, atom_type]];  % Append platelet points to coordinate matrix
    end
end

% ==========================================================
% Section 6: Create 3D Mesh (400x400x400)
% ==========================================================
mesh_space = zeros(grid_size, grid_size, grid_size); % Initialize the 3D mesh
scaling_factor = max(max(abs(Points))); % Scaling factor to fit into the 400x400x400 grid

% Populate the 3D mesh with Fibrin and Platelet atoms
for i = 1:size(coordinate_matrix, 1)
    x = coordinate_matrix(i, 1);
    y = coordinate_matrix(i, 2);
    z = coordinate_matrix(i, 3);
    atom_type = coordinate_matrix(i, 4);

    % Convert coordinates to indices in the 3D matrix
    x_idx = round((x / scaling_factor) * (grid_size / 2)) + grid_size / 2;
    y_idx = round((y / scaling_factor) * (grid_size / 2)) + grid_size / 2;
    z_idx = round((z / scaling_factor) * (grid_size / 2)) + grid_size / 2;

    % Ensure indices are within bounds
    x_idx = max(min(x_idx, grid_size), 1);
    y_idx = max(min(y_idx, grid_size), 1);
    z_idx = max(min(z_idx, grid_size), 1);

    % Set the value in the mesh: 1 for Fibrin, 2 for Platelet
    if atom_type == 1
        mesh_space(x_idx, y_idx, z_idx) = 1; % Fibrin
    elseif atom_type == 2
        mesh_space(x_idx, y_idx, z_idx) = 2; % Platelet
    end
end

% ==========================================================
% Section 7: Plot Final Network with Bezier Curves and Platelets
% ==========================================================
% ==========================================================
% Section 7: Plot Final Network with Bezier Curves, Platelets, and Sphere Inclusions
% ==========================================================
figure;
hold on;

% Plot Bezier curves for bonds
for i = 1:size(bonds, 1)
    P0 = Points(bonds(i, 1), :);  % Start point of the bond
    P2 = Points(bonds(i, 2), :);  % End point of the bond
    midpoint = (P0 + P2) / 2;  % Midpoint for control point
    random_curvature = (rand(1, 3) - 0.5) * 3;  % Random curvature
    P1 = midpoint + random_curvature;  % Control point for Bezier curve

    t = linspace(0, 1, num_bond_points);  % Parameter for Bezier curve
    bezier_curve = (1 - t').^2 * P0 + 2 * (1 - t') .* t' * P1 + t'.^2 * P2;  % Bezier curve formula

    % Plot the Bezier curve
    plot3(bezier_curve(:,1), bezier_curve(:,2), bezier_curve(:,3), 'r-', 'LineWidth', 1);
end

% Plot platelets (spheres)
for i = 1:length(platelet_centers)
    [x_sphere, y_sphere, z_sphere] = sphere(20);  % Create a sphere
    x_sphere = x_sphere * (platelet_diameter / 2) + platelet_centers(i, 1);
    y_sphere = y_sphere * (platelet_diameter / 2) + platelet_centers(i, 2);
    z_sphere = z_sphere * (platelet_diameter / 2) + platelet_centers(i, 3);
    surf(x_sphere, y_sphere, z_sphere, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'g');
end

% Plot the sphere inclusions
for i = 1:length(sphere_radii)
    [x_sphere, y_sphere, z_sphere] = sphere(20);  % Create a sphere with resolution 20
    x_sphere = x_sphere * sphere_radii(i) + sphere_centers(i, 1);
    y_sphere = y_sphere * sphere_radii(i) + sphere_centers(i, 2);
    z_sphere = z_sphere * sphere_radii(i) + sphere_centers(i, 3);
    surf(x_sphere, y_sphere, z_sphere, 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'FaceColor', 'b');
end

axis equal;
title('Final Network with Bezier Curves, Platelets, and Sphere Inclusions');
hold off;
view(3);


% ==========================================================
% Section 8: Save Results
% ==========================================================
save('clot_data.mat', 'coordinate_matrix', 'mesh_space');
disp('Data saved to clot_data.mat.');

% ==========================================================
% Section 9: Helper Functions
% ==========================================================
function [bonds, bond_lengths] = generate_initial_bonds(Points)
% Generate initial bonds with random number of connections (3 to 4) per node
numPoints = size(Points, 1);
bonds = [];
connection_counts = zeros(numPoints, 1); % Track the number of connections per point
target_connections = randi([3, 4], numPoints, 1); % Randomly assign 3 or 4 connections to each point

for i = 1:numPoints
    distances = sqrt(sum((Points - Points(i, :)).^2, 2));
    [~, idx] = sort(distances);
    nearest_neighbors = idx(2:end); % Nearest neighbors, excluding self

    for j = 1:length(nearest_neighbors)
        if connection_counts(i) < target_connections(i) && connection_counts(nearest_neighbors(j)) < target_connections(nearest_neighbors(j))
            bonds = [bonds; i, nearest_neighbors(j)];
            connection_counts(i) = connection_counts(i) + 1;
            connection_counts(nearest_neighbors(j)) = connection_counts(nearest_neighbors(j)) + 1;
        end

        % Stop adding connections if the point has reached its target connections
        if connection_counts(i) >= target_connections(i)
            break;
        end
    end
end
bond_lengths = compute_bond_lengths(bonds, Points);
end

function bond_lengths = compute_bond_lengths(bonds, Points)
% Compute the lengths of bonds between points
bond_lengths = zeros(size(bonds, 1), 1);
for i = 1:size(bonds, 1)
    bond_lengths(i) = norm(Points(bonds(i, 1), :) - Points(bonds(i, 2), :));
end
end

function forces = compute_forces(bonds, bond_lengths, d0_fibrin, k_fibrin, Points)
% Compute forces for each point based on bond lengths
forces = zeros(size(Points));
for i = 1:size(bonds, 1)
    bond_vec = Points(bonds(i, 2), :) - Points(bonds(i, 1), :);
    force_magnitude = k_fibrin * (bond_lengths(i) - d0_fibrin);
    unit_vec = bond_vec / norm(bond_vec);
    forces(bonds(i, 1), :) = forces(bonds(i, 1), :) + force_magnitude * unit_vec;
    forces(bonds(i, 2), :) = forces(bonds(i, 2), :) - force_magnitude * unit_vec;
end
end

function Points = move_points(Points, forces, step_size)
% Move points based on computed forces
Points = Points + step_size * forces;
end
