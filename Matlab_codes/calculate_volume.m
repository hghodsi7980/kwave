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