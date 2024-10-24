% ==========================================================
% Fibrin Network Optimization for Mechanical Stability
% ==========================================================
clear; clc;

% -------------------
% Parameters
% -------------------
numPoints = 5000; % Number of points in the fibrin network
max_iterations = 1000; % Max number of optimization iterations
energy_threshold = 1e-3; % Energy threshold for optimization stop
k_fibrin = 1; % Stiffness constant for fibrin fibers (in nN/µm)
step_size = 0.01; % Step size for moving points during optimization

% ==========================================================
% Section 1: Initialization of Fibrin Network
% ==========================================================
% Generate random points within an ellipsoid as fibrin network points
a = round(max(20, abs(normrnd(75, 25))));
b = round(max(20, abs(normrnd(75, 25))));
c = round(max(20, abs(normrnd(75, 25))));

% Calculate the volume of the ellipsoid
ellipsoid_volume = (4/3) * pi * a * b * c;

% Calculate the equilibrium bond length based on volume and number of points
d0_fibrin = (ellipsoid_volume / numPoints)^(1/3);

disp(['Calculated equilibrium bond length (d0_fibrin): ', num2str(d0_fibrin)]);

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

% Generate initial bonds between points (Delaunay triangulation or nearest neighbors)
[bonds, bond_lengths] = generate_initial_bonds(Points);

% ==========================================================
% Section 2: Energy Function and Optimization Loop
% ==========================================================
compute_energy = @(distances) sum(k_fibrin * (distances - d0_fibrin).^2);

% Optimization loop
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
    
    % Step 2: Vary bond configurations to reduce energy
    for i = 1:size(bonds, 1)
        % Pass Points to vary_bonds function
        new_bonds = vary_bonds(bonds, Points); 
        new_bond_lengths = compute_bond_lengths(new_bonds, Points);
        new_energy = compute_energy(new_bond_lengths);
        
        % Accept the new configuration if energy is reduced
        if new_energy < current_energy
            bonds = new_bonds;
            bond_lengths = new_bond_lengths;
            current_energy = new_energy;
        end
    end
    
    % Step 3: Move points to reduce discrepancies in bond lengths
    forces = compute_forces(bonds, bond_lengths, d0_fibrin, k_fibrin, Points); % Calculate forces based on bond lengths
    Points = move_points(Points, forces, step_size); % Move points based on forces
    
    % Recompute bond lengths after moving points
    bond_lengths = compute_bond_lengths(bonds, Points);
    
    % Display current energy for debugging
    disp(['Iteration ', num2str(iteration), ', Energy: ', num2str(current_energy)]);
end

% ==========================================================
% Section 3: Visualization
% ==========================================================
% Plot the optimized fibrin network
figure;
scatter3(Points(:,1), Points(:,2), Points(:,3), 'filled');
hold on;
for i = 1:size(bonds, 1)
    plot3([Points(bonds(i, 1), 1), Points(bonds(i, 2), 1)], ...
          [Points(bonds(i, 1), 2), Points(bonds(i, 2), 2)], ...
          [Points(bonds(i, 1), 3), Points(bonds(i, 2), 3)], 'r-');
end
axis equal;
title('Optimized Fibrin Network');
hold off;

% ==========================================================
% Section 4: Helper Functions
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

function new_bonds = vary_bonds(bonds, Points)
    % Vary bond configurations (e.g., swap or reconnect bonds)
    new_bonds = bonds;
    % Randomly vary the connections
    if rand() < 0.5
        swap_idx = randi([1 size(bonds, 1)], 1, 2);
        new_bonds([swap_idx(1), swap_idx(2)], :) = new_bonds([swap_idx(2), swap_idx(1)], :);
    else
        % Optionally reconnect one bond randomly
        reconnect_idx = randi([1 size(bonds, 1)], 1);
        new_bonds(reconnect_idx, 2) = randi([1 size(Points, 1)], 1);
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
