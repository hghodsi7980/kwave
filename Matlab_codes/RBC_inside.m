% Step 1: Generate the sphere with radius = 50
radius_sphere = 50;
[x_sphere, y_sphere, z_sphere] = sphere(100);  % Higher resolution sphere
x_sphere = x_sphere * radius_sphere;
y_sphere = y_sphere * radius_sphere;
z_sphere = z_sphere * radius_sphere;

% Plot the sphere
figure;
surf(x_sphere, y_sphere, z_sphere, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
axis equal;
hold on;

% Step 2: Define RBC parameters
scaling_factor_RBC = 1;
num_RBCs = 700;  % Increase number of RBCs to test
rbc_diameter = 8 / scaling_factor_RBC;  % Approximate RBC size

% Step 3: Randomly generate positions inside the sphere and allow adjacent placement
rbc_positions = [];
min_distance = rbc_diameter;  % Distance between RBCs (just beside, no intersection)

% Place the first RBC randomly inside the sphere
theta = 2 * pi * rand();
phi = acos(2 * rand() - 1);
r = radius_sphere * rand()^(1/3);  % Cubic root for uniform distribution inside sphere
rbc_center = [r * sin(phi) * cos(theta), r * sin(phi) * sin(theta), r * cos(phi)];
rbc_positions = [rbc_positions; rbc_center];

% Step 4: Place remaining RBCs adjacent to existing ones
for i = 2:num_RBCs
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
        if norm(new_rbc_center) <= radius_sphere
            % Ensure the new RBC doesn't intersect with any existing RBCs
            distances = sqrt(sum((rbc_positions - new_rbc_center).^2, 2));
            if all(distances > min_distance)
                rbc_positions = [rbc_positions; new_rbc_center];  % Add new RBC position
                break;
            end
        end
    end
end

% Step 5: Calculate RBC volumes using convex hulls
total_rbc_volume = 0;  % Initialize total RBC volume
for i = 1:num_RBCs
    center_point_t = rbc_positions(i, :);
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
        trisurf(K, rbc_points(:, 1), rbc_points(:, 2), rbc_points(:, 3), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'r');  % Plot RBC surface
        
        % Calculate and accumulate volume of the convex hull
        [~, rbc_volume] = convhull(rbc_points);
        total_rbc_volume = total_rbc_volume + rbc_volume;
    end
end

% Step 6: Calculate sphere volume
sphere_volume = (4/3) * pi * (radius_sphere^3);

% Step 7: Calculate the filling factor
filling_factor = total_rbc_volume / sphere_volume;

disp(['Total RBC volume: ', num2str(total_rbc_volume)]);
disp(['Sphere volume: ', num2str(sphere_volume)]);
disp(['Filling factor: ', num2str(filling_factor)]);

hold off;
title([num2str(num_RBCs), ' RBCs in Random Positions Inside a Sphere with Convex Hulls']);
axis equal;
grid on;
