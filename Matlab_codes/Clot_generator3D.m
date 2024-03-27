clear
close all
hold on
grid on
grid minor
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

% trisurf(C,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3), ...
% 'EdgeColor','black','FaceColor','none')

% Set of coordinates to check
number_of_overall_points = round(max([a,b,c])*2.5*abs((2+randn())));
concentration_factor = 0.25*abs((2+randn()));
% Calculate the center point of the surface
center_point = [mean(DT.Points(C, 1)), mean(DT.Points(C, 2)), mean(DT.Points(C, 3))];
x_center = center_point(1);
y_center = center_point(2);
z_center = center_point(3);

x_gaussian = normrnd(x_center, concentration_factor*a, number_of_overall_points, 1); % Adjust standard deviation as needed
y_gaussian = normrnd(y_center, concentration_factor*b, number_of_overall_points, 1);
z_gaussian = normrnd(z_center, concentration_factor*c, number_of_overall_points, 1);

coordinates_to_check = [x_gaussian,y_gaussian,z_gaussian];

% Check if points are inside or on the surface
in_or_on_ellipsoid = ~isnan(pointLocation(DT, coordinates_to_check));

% Choose the inside points
inside_points = coordinates_to_check(in_or_on_ellipsoid, :);

% Create Delaunay triangulation of points inside the ellipsoid
DT_inside = delaunayTriangulation(inside_points);


% Extract tetrahedral elements from Delaunay triangulation
tetrahedra = DT_inside.ConnectivityList;

% Initialize variables for average bond length calculation
total_bond_length = 0;
num_bonds = 0;

% Initialize matrix to store all coordinates and material indices
all_coordinates = [];

% Generate waitbar
holder = waitbar(0, 'Progress RBC + Fibrin'); % Create a progress bar window
% Number_of bonds per point
bond_per_point = zeros(size(inside_points,1),1);
% Loop through each tetrahedron
for i = 1:size(tetrahedra, 1)
    v = tetrahedra(i, :);
    % Update waitbar
    waitbar(i / size(tetrahedra, 1), holder, sprintf('Progress RBC + Fibrin: %d%%',round((i / size(tetrahedra, 1)) * 100)));
    % Compute centroid (center point) of the tetrahedron
    center_point = mean(inside_points(v, :), 1);

    % Plot center point
    distance_to_origin = sqrt(sum(center_point.^2));
    show_RBC = abs((distance_to_origin/concentration_factor)*randn()) <= max([a, b, c]);

    if show_RBC
        % Define constants for RBCs
        d = 8; % in meters
        b = 1; % in meters
        h = 2.12; % in meters
        % Calculate P, Q, and R for RBCs
        P = -(d^2/2) + (h^2/2) * ((d^2/b^2) - 1) - h^2/2 * ((d^2/b^2) - 1) * sqrt(1 - (b^2/h^2));
        Q = P * (d^2/b^2) + (b^2/4) * (d^4/b^4 - 1);
        R = -P * (d^2/4) - d^4/16;
        [x_rbc,y_rbc,z_rbc] = meshgrid(-10+center_point(1):0.5:10+ ...
            center_point(1),-10+center_point(2):0.5:10+center_point(2),...
            -10+center_point(3):0.5:10+center_point(3));
        cell_rotations = rand(1, 3) *pi;
        % Compute transformed coordinates for RBCs
        x_rbc_rot = x_rbc-center_point(1);
        y_rbc_rot = y_rbc-center_point(2);
        z_rbc_rot = z_rbc-center_point(3);
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

    % Calculate bond length
    bond_length = sqrt(sum((inside_points(v(1), :) - inside_points(v(2), :)).^2));
    % Check if current bond length is shorter or equal to average
    %check the current number of bond per points
    % Calculate X, Y, and Z using inside_points
    edges = [
        v(1), v(2);
        v(1), v(3);
        v(1), v(4);
        v(2), v(3);
        v(2), v(4);
        v(3), v(4)
        ];
    line_coordinates = []; % Initialize line coordinates
    X = inside_points(edges(:, 1), 1);
    Y = inside_points(edges(:, 1), 2);
    Z = inside_points(edges(:, 1), 3);
    pt1 = [X(1), Y(1), Z(1)];
    pt2 = [X(4), Y(4), Z(4)];
    pt3 = [X(6), Y(6), Z(6)];
    point_index = find(inside_points == pt2);
    if num_bonds > 0 && bond_length/2 <= (total_bond_length / num_bonds) && bond_per_point(point_index(1)) < 4
        plot3(X',Y',Z','b');
        bond_per_point(point_index(1)) = bond_per_point(point_index(1))+2;
        % Linearly interpolate points between pt1 and pt2
        interp_pts = [];
        for t = linspace(0, 1, 100)
            interp_pt = (1 - t) * pt1 + t * pt2;
            interp_pts = [interp_pts; interp_pt];
        end
        % Append interpolated points to line_coordinates
        line_coordinates = [line_coordinates; interp_pts];

        % Linearly interpolate points between pt2 and pt3
        interp_pts = [];
        for t = linspace(0, 1, 50)
            interp_pt = (1 - t) * pt2 + t * pt3;
            interp_pts = [interp_pts; interp_pt];
        end
        % Append interpolated points to line_coordinates
        line_coordinates = [line_coordinates; interp_pts];

        % Generate material indices for each point
        material_index_line = ones(size(line_coordinates, 1), 1) * 2; % Material index for lines
        % Appending to line coordinates
        line_coordinates = [line_coordinates, material_index_line];
        % Append line coordinates to all_coordinates
        all_coordinates = [all_coordinates; line_coordinates];
    end


    % Increment total bond length and number of bonds
    total_bond_length = total_bond_length + bond_length;
    num_bonds = num_bonds + 1;
end

close(holder);

% Randomly exclude some points
exclude_percentage_platelet = rand()*0.5;
num_points_to_exclude = round(exclude_percentage_platelet * size(inside_points, 1));
indices_to_exclude = randperm(size(inside_points, 1), num_points_to_exclude);
inside_points(indices_to_exclude, :) = [];



platelet_radius = 1.5;
for i = 1:size(inside_points,1)

    % Generate surface coordinates for the sphere
    [sx, sy, sz] = sphere;
    sx = sx * platelet_radius + inside_points(i, 1);
    sy = sy * platelet_radius + inside_points(i, 2);
    sz = sz * platelet_radius + inside_points(i, 3);
    surf(sx,sy,sz,'FaceColor', 'green', 'EdgeColor', 'none');
    % Append surface sphere coordinates to all_coordinates
    sphere_surface_coordinates = [sx(:), sy(:), sz(:), ones(size(sx(:)))*3]; % Material index for spheres
    all_coordinates = [all_coordinates; sphere_surface_coordinates]; %#ok<*AGROW>

    % Generate inside coordinates for the sphere
    num_inside_points = 500; % Adjust as needed
    sphere_center = inside_points(i, :);
    sphere_radius = platelet_radius;

    % Generate points uniformly distributed inside the sphere
    rand_points = rand(num_inside_points, 3) - 0.5; % Points in [-0.5, 0.5]
    rand_points = rand_points ./ sqrt(sum(rand_points.^2, 2)); % Normalize to unit vectors
    inside_sphere_coordinates = sphere_center + sphere_radius * rand_points;

    % Append inside sphere coordinates to all_coordinates
    sphere_inside_coordinates = [inside_sphere_coordinates, ones(size(inside_sphere_coordinates, 1), 1)*3]; % Material index for spheres
    all_coordinates = [all_coordinates; sphere_inside_coordinates];
end

index = 3;

% Sort the coordinates matrix based on first x, then y, and finally z
all_coordinates_sorted = sortrows(single(all_coordinates), [4,1,2,3]);
% % Save sorted coordinates to a text file
writematrix(all_coordinates_sorted, sprintf('all_coordinates_sorted_%d.txt',index), 'Delimiter', 'tab');

coordinates = all_coordinates_sorted(:,1:3);
materials = all_coordinates_sorted(:,4);


% Define the size of the mesh space
mesh_size = 2000;

% Create an empty mesh space
mesh_space = zeros(mesh_size, mesh_size, mesh_size, 'single');

% Translate coordinates to positive quadrant
coordinates = coordinates - min(coordinates) + 1;

% Scale coordinates to fit within the mesh size
scaling_factor = mesh_size / max(max(coordinates));
coordinates_scaled = coordinates * scaling_factor;

% Ensure all coordinates are integers
coordinates_scaled = round(coordinates_scaled);

% Assign material indices to the mesh space
for i = 1:size(coordinates_scaled, 1)
    x = coordinates_scaled(i, 1);
    y = coordinates_scaled(i, 2);
    z = coordinates_scaled(i, 3);

    % Check if the coordinates are within the mesh space
    if x >= 1 && x <= mesh_size && y >= 1 && y <= mesh_size && z >= 1 && z <= mesh_size
        % Assign the material index to the corresponding voxel
        mesh_space(x, y, z) = materials(i);
    end
end

% Save the mesh space to a file
save('mesh_space.mat', 'mesh_space', '-v7.3');

