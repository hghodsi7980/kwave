% ==========================================================
% Fibrin Network with Strain Energy Relaxation, Bezier Curves, and Random Non-Intersecting Spheres
% ==========================================================
clear; clc; close all

% -------------------
% Parameters
% -------------------
numPoints = 5000; % Number of points in the fibrin network
max_iterations = 1000; % Max number of optimization iterations
energy_threshold = 1e-3; % Energy threshold for optimization stop
k_fibrin = 1; % Stiffness constant for fibrin fibers (in nN/µm)
step_size = 0.01; % Step size for moving points during optimization
numSpheres = randi(50); % Number of random spheres (inclusions)
max_diameter_factor = 0.5; % Maximum sphere diameter as a fraction of the largest dimension

% ==========================================================
% Section 1: Initialization of Fibrin Network
% ==========================================================
% Generate random points within an ellipsoid as fibrin network points
a = round(max(20, abs(normrnd(75, 25))));
b = round(max(20, abs(normrnd(75, 25))));
c = round(max(20, abs(normrnd(75, 25))));


% Generate random points on the ellipsoid
theta = 2 * pi * rand(numPoints, 1);
phi = pi * rand(numPoints, 1);
x_ellipsoid = a * sin(phi) .* cos(theta);
y_ellipsoid = b * sin(phi) .* sin(theta);
z_ellipsoid = c * cos(phi);

% Generate random displacements
x_displacement = a * rand(numPoints, 1);
y_displacement = b * rand(numPoints, 1);
z_displacement = c * rand(numPoints, 1);
Points = [x_ellipsoid + x_displacement, y_ellipsoid + y_displacement, z_ellipsoid + z_displacement];

% Generate initial bonds between points (Nearest neighbors)
[bonds, bond_lengths] = generate_initial_bonds(Points);

% ==========================================================
% Section 2: Strain Energy Relaxation of the Network
% ==========================================================
% Calculate the volume of the ellipsoid

DT = delaunayTriangulation(Points);
[C, ellipsoid_volume] = convexHull(DT);

% Calculate the equilibrium bond length based on volume and number of points
d0_fibrin = (ellipsoid_volume / numPoints)^(1/3);

disp(['Calculated equilibrium bond length (d0_fibrin): ', num2str(d0_fibrin)]);

compute_energy = @(distances) sum(k_fibrin * (distances - d0_fibrin).^2);

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
end

% ==========================================================
% Plot the Relaxed Network with Bezier Curves
% ==========================================================
figure;
scatter3(Points(:,1), Points(:,2), Points(:,3), 'filled');
hold on;
% Plot bonds as Bezier curves with random curvature
for i = 1:size(bonds, 1)
    P0 = Points(bonds(i, 1), :);  % Start point of the bond
    P2 = Points(bonds(i, 2), :);  % End point of the bond

    % Midpoint of the bond (for generating the control point)
    midpoint = (P0 + P2) / 2;

    % Generate a random curvature by perturbing the midpoint
    if rand() < 0.9  % 90% chance for small curvature, 10% for large curvature
        random_curvature = (rand(1, 3) - 0.5) * 3;  % Small curvature
    else
        random_curvature = (rand(1, 3) - 0.5) * 15;  % Large curvature
    end
    P1 = midpoint + random_curvature;  % Control point

    % Generate the Bezier curve
    t = linspace(0, 1, 50);  % Parameter t for Bezier curve
    bezier_curve = (1 - t).^2' * P0 + 2 * (1 - t)' .* t' * P1 + t.^2' * P2;

    % Plot the Bezier curve
    plot3(bezier_curve(:,1), bezier_curve(:,2), bezier_curve(:,3), 'r-', 'LineWidth', 1);
end
axis equal;
title('Relaxed Fibrin Network with Bezier Curves');
hold off;

% ==========================================================
% Section 3: Add Random Spheres (Non-Intersecting in Extended Network Volume)
% ==========================================================
max_dimension = max([a, b, c]);  % Largest dimension of the initial shape
max_diameter = max_dimension * max_diameter_factor;

% Calculate the centroid of the fibrin network (center of the network)
centroid = mean(Points);

% Standard deviation for the normal distribution of sphere centers
std_dev_factor = 0.5; % Adjust this factor to control how clustered the spheres are around the center

sphere_centers = [];
sphere_radii = [];

for i = 1:numSpheres
    while true
        % Generate random sphere center using normal distribution around the network centroid
        center = centroid + std_dev_factor * [std(Points(:,1)) * randn(), ...
                                              std(Points(:,2)) * randn(), ...
                                              std(Points(:,3)) * randn()];
        radius = max(max_diameter / 6 , rand() * max_diameter / 2);

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


% ==========================================================
% Plot the Network with Spheres and Bezier Curves
% ==========================================================
figure;
scatter3(Points(:,1), Points(:,2), Points(:,3), 'filled');
hold on;
% Plot Bezier curves for bonds
for i = 1:size(bonds, 1)
    P0 = Points(bonds(i, 1), :);  % Start point of the bond
    P2 = Points(bonds(i, 2), :);  % End point of the bond

    % Midpoint of the bond (for generating the control point)
    midpoint = (P0 + P2) / 2;

    % Generate a random curvature by perturbing the midpoint
    if rand() < 0.9  % 90% chance for small curvature, 10% for large curvature
        random_curvature = (rand(1, 3) - 0.5) * 3;  % Small curvature
    else
        random_curvature = (rand(1, 3) - 0.5) * 15;  % Large curvature
    end
    P1 = midpoint + random_curvature;  % Control point

    % Generate the Bezier curve
    t = linspace(0, 1, 50);  % Parameter t for Bezier curve
    bezier_curve = (1 - t).^2' * P0 + 2 * (1 - t)' .* t' * P1 + t.^2' * P2;

    % Plot the Bezier curve
    plot3(bezier_curve(:,1), bezier_curve(:,2), bezier_curve(:,3), 'r-', 'LineWidth', 1);
end

% Plot spheres
for i = 1:length(sphere_radii)
    [x_sphere, y_sphere, z_sphere] = sphere(20);  % Create a sphere
    x_sphere = x_sphere * sphere_radii(i) + sphere_centers(i, 1);
    y_sphere = y_sphere * sphere_radii(i) + sphere_centers(i, 2);
    z_sphere = z_sphere * sphere_radii(i) + sphere_centers(i, 3);
    surf(x_sphere, y_sphere, z_sphere, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
axis equal;
title('Fibrin Network with Random Spheres and Bezier Curves');
hold off;

% ==========================================================
% Section 4: Remove Nodes and Fibers Inside the Spheres
% ==========================================================
inside_sphere = false(size(Points, 1), 1);  % To mark nodes inside spheres

% Check each point if it lies inside any of the spheres
for i = 1:size(Points, 1)
    for j = 1:length(sphere_radii)
        dist_to_sphere_center = norm(Points(i, :) - sphere_centers(j, :));
        if dist_to_sphere_center < sphere_radii(j)
            inside_sphere(i) = true;
            break;
        end
    end
end

% Remove nodes that are inside spheres
Points(inside_sphere, :) = [];

% Update bond indices after removing nodes
new_indices = find(~inside_sphere);  % Indices of the remaining nodes
mapping = zeros(size(inside_sphere));  % Create a mapping from old to new indices
mapping(new_indices) = 1:length(new_indices);  % Fill in the mapping

% Update bonds based on the new mapping
bonds = bonds(~ismember(bonds(:,1), find(inside_sphere)) & ~ismember(bonds(:,2), find(inside_sphere)), :);
bonds = [mapping(bonds(:, 1)), mapping(bonds(:, 2))];  % Re-map bond indices to the new Points array

% ==========================================================
% Section 5: Squeezing the Network Based on Sphere Influence
% ==========================================================
for i = 1:size(Points, 1)
    for j = 1:length(sphere_radii)
        % Calculate distance from point to sphere center
        dist_to_center = norm(Points(i, :) - sphere_centers(j, :));
        if dist_to_center < sphere_radii(j) * 2  % Influence zone
            % Move point outward away from the sphere center
            direction = (Points(i, :) - sphere_centers(j, :)) / dist_to_center;
            influence_factor = (sphere_radii(j) * 2 - dist_to_center) / (sphere_radii(j) * 2);
            Points(i, :) = Points(i, :) + influence_factor * direction * 5;  % Move the point outward
        end
    end
end

% ==========================================================
% Plot the Squeezed Network with Bezier Curves
% ==========================================================
figure;
scatter3(Points(:,1), Points(:,2), Points(:,3), 'filled','g');
hold on;
% Plot Bezier curves for bonds
for i = 1:size(bonds, 1)
    P0 = Points(bonds(i, 1), :);  % Start point of the bond
    P2 = Points(bonds(i, 2), :);  % End point of the bond

    % Midpoint of the bond (for generating the control point)
    midpoint = (P0 + P2) / 2;

    % Generate a random curvature by perturbing the midpoint
    if rand() < 0.9  % 90% chance for small curvature, 10% for large curvature
        random_curvature = (rand(1, 3) - 0.5) * 3;  % Small curvature
    else
        random_curvature = (rand(1, 3) - 0.5) * 15;  % Large curvature
    end
    P1 = midpoint + random_curvature;  % Control point

    % Generate the Bezier curve
    t = linspace(0, 1, 50);  % Parameter t for Bezier curve
    bezier_curve = (1 - t).^2' * P0 + 2 * (1 - t)' .* t' * P1 + t.^2' * P2;

    % Plot the Bezier curve
    plot3(bezier_curve(:,1), bezier_curve(:,2), bezier_curve(:,3), 'r-', 'LineWidth', 1);
end

% Plot spheres
for i = 1:length(sphere_radii)
    [x_sphere, y_sphere, z_sphere] = sphere(20);  % Recreate spheres
    x_sphere = x_sphere * sphere_radii(i) + sphere_centers(i, 1);
    y_sphere = y_sphere * sphere_radii(i) + sphere_centers(i, 2);
    z_sphere = z_sphere * sphere_radii(i) + sphere_centers(i, 3);
    surf(x_sphere, y_sphere, z_sphere, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
axis equal;
title('Squeezed Fibrin Network with Bezier Curves and Spheres');
hold off;

% ==========================================================
% Section 6: Helper Functions
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
