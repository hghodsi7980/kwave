% Define control points for a closed 3D B-spline volume
control_points = zeros(3, 3, 3, 3); % 3x3x3x3 grid of control points

% Define control points for the first layer
control_points(:, :, 1, 1) = [0 100 200; 0 100 200; 0 100 200];

% Copy control points for the remaining layers with slight modifications
for k = 2:3
    control_points(:, :, :, k) = control_points(:, :, :, 1) + (k - 1) * 100;
end

% Create knots for the B-spline volume
knots = [0 0 0 1 1 1]; % Example knots, adjust as needed

% Create B-spline surfaces in each direction
spline_x = spmak(knots, control_points(1, :, :, :));
spline_y = spmak(knots, control_points(2, :, :, :));
spline_z = spmak(knots, control_points(3, :, :, :));

% Evaluate the B-spline surfaces at a grid of points to create the volume
n_eval_points = 50;
eval_points = linspace(0, 1, n_eval_points);
[X, Y, Z] = ndgrid(eval_points, eval_points, eval_points);
points = [X(:), Y(:), Z(:)];

% Evaluate B-spline surfaces at the grid of points
evaluated_x = fnval(spline_x, points(:, 1));
evaluated_y = fnval(spline_y, points(:, 2));
evaluated_z = fnval(spline_z, points(:, 3));

% Reshape evaluated points to match the grid
evaluated_x_grid = reshape(evaluated_x, size(X));
evaluated_y_grid = reshape(evaluated_y, size(Y));
evaluated_z_grid = reshape(evaluated_z, size(Z));

% Visualize the evaluated B-spline volume
scatter3(evaluated_x_grid(:), evaluated_y_grid(:), evaluated_z_grid(:), 'b', '.');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;
grid on;
