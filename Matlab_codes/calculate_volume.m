function total_volume = calculate_volume(matrix)
% Convert the matrix to a binary matrix where 1 represents material
binary_matrix = matrix > 0;

% Define a smaller structuring element (allowing for less than 6-connectivity)
% A spherical or cubic element to allow more lenient filling
se = strel('sphere', 1);  % This creates a spherical structuring element with radius 1
% Alternatively, you can use a cubic structuring element like:
% se = strel('cube', 2);  % This allows a more relaxed connection

% Perform dilation to fill in gaps between materials with the custom structuring element
filled_matrix = imdilate(binary_matrix, se);

% Convert the filled binary matrix back to the original matrix values
filled_matrix_values = filled_matrix .* matrix;

% Calculate the total volume of non-zero elements (sum of all 1's in the filled matrix)
total_volume = sum(filled_matrix(:));

end
% mesh_space(mesh_space == 0) = 4;
% [x, y, z] = ind2sub(size(mesh_space), find(mesh_space > 0));
% values = mesh_space(mesh_space > 0);
% colors = zeros(length(values), 3); % initialize color array
% colors(values == 4, :) = repmat([1 0 0], sum(values == 4), 1); % Red
% colors(values == 1, :) = repmat([1 0 0], sum(values == 1), 1); % Red
% colors(values == 2, :) = repmat([0 1 0], sum(values == 2), 1); % Green
% colors(values == 3, :) = repmat([0 0 1], sum(values == 3), 1); % Blue
% scatter3(x, y, z, 1, colors, 'filled');
% 
% grid on;
% axis equal
% view(125,25);