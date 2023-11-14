% fill an area specified by a mask with RBC cells in an Nx * Ny matrix
% the output is a mask in wich all of the RBCs areas are 1 while other
% indices are zero
% RBC_radius is a number related to Nx and Ny
% licenced to Hamed Ghodsi
% first version 16-10-2023
function output_mask = fill_with_RBC(Nx, Ny, filling_mask, RBC_count, RBC_radius)
    output_mask = zeros(Nx, Ny);
    occupied_regions = zeros(Nx, Ny);

    [center_x, center_y] = find(filling_mask);
    available_coords = [center_x, center_y];
    
    for k = 1:RBC_count
%         if isempty(available_coords)
%             break; % No more available coordinates
%         end
        disp(k*100/RBC_count)
        % Randomly select an available coordinate
        random_index = randi(size(available_coords, 1));
        center_x = available_coords(random_index, 1);
        center_y = available_coords(random_index, 2);

        % Remove the selected coordinate from the list
        available_coords(random_index, :) = [];

        % Set the center element to 1
        output_mask(center_x, center_y) = 1;

        % Create the circular region with radius "RBC_radius" without overlap
        for i = 1:Nx
            for j = 1:Ny
                if (i - center_x)^2 + (j - center_y)^2 <= RBC_radius^2
                    output_mask(i, j) = 1;
                    occupied_regions(i, j) = 0; % Mark the region as occupied
                end
            end
        end
    end
end
