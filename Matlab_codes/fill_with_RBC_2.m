function output_mask = fill_with_RBC_2(Nx, Ny, filling_mask, RBC_count, RBC_radius)
output_mask = zeros(Nx, Ny);
occupied_regions = zeros(Nx, Ny);

[center_x, center_y] = find(filling_mask);
available_coords = [center_x, center_y];

% Calculate the number of rows and columns in the grid
num_rows = ceil(sqrt(RBC_count));
num_cols = ceil(RBC_count / num_rows);

% Calculate the grid spacing
grid_spacing_x = floor(Nx / num_cols);
grid_spacing_y = floor(Ny / num_rows);

for row = 1:num_rows
    for col = 1:num_cols
        % Check if there are available coordinates
        if isempty(available_coords)
            break; % No more available coordinates
        end

        % Get the next available coordinate
        coord = available_coords(1, :);

        % Remove the selected coordinate from the list
        available_coords(1, :) = [];

        % Calculate the center of the circle based on the grid
        center_x = (row - 1) * grid_spacing_x + floor(grid_spacing_x / 2);
        center_y = (col - 1) * grid_spacing_y + floor(grid_spacing_y / 2);

        % Adjust the center based on the available coordinate
        center_x = center_x + coord(1) - floor(Nx / 2);
        center_y = center_y + coord(2) - floor(Ny / 2);

        % Set the center element to 1
        output_mask(coord(1), coord(2)) = 1;

        % Create the circular region with radius "RBC_radius" without overlap
        for i = 1:Nx
            for j = 1:Ny
                if (i - coord(1))^2 + (j - coord(2))^2 <= RBC_radius^2
                    output_mask(i, j) = 1;
                    occupied_regions(i, j) = 1; % Mark the region as occupied
                end
            end
        end
    end
end
end
