% fill an area specified by a mask with Fibrin strands in an Nx * Ny matrix
% the output is a mask in wich all of the Fibrin areas are 2 while other
% indices are zero
% licenced to Hamed Ghodsi
% first version 1-11-2023
function output_mask = fill_with_fibrin(Nx, Ny, filling_mask, Fibrin_count, max_fibrin_length)
output_mask = zeros(Nx, Ny);
h = waitbar(0, 'Progress'); % Create a progress bar window
for k = 1:Fibrin_count

    waitbar(k / Fibrin_count, h, sprintf('Fill with Fibrin Progress: %d%%', round(k*100/Fibrin_count)));
    % pause(0.1); % Optional: Add a delay to slow down the progress for demonstration purposes
    % Randomly select a position within the filling_mask
    [x, y] = find(filling_mask);
    random_index = randi(length(x));
    start_x = x(random_index)+randi([-95 95]);
    start_y = y(random_index)+randi([-95 95]);

    % Randomly generate a length between 1 and max_fibrin_length
    fibrin_length = max(1, round(rand() * max_fibrin_length));

    % Generate a random angle between 0 and 90 degrees
    angle = (pi / 2) * rand();

    % Calculate the endpoints of the line
    end_x = start_x + round(fibrin_length * cos(angle));
    end_y = start_y + round(fibrin_length * sin(angle));

    % Ensure that the line stays within the boundaries
    % start_x = max(1, min(Nx, start_x));
    % start_y = max(1, min(Ny, start_y));
    % end_x = max(1, min(Nx, end_x));
    % end_y = max(1, min(Ny, end_y));

    % Create the line by drawing pixels between start and end points
    x_values = linspace(start_x, end_x, fibrin_length);
    y_values = linspace(start_y, end_y, fibrin_length);

    for i = 1:fibrin_length
        x = round(x_values(i));
        y = round(y_values(i));
        % if filling_mask(x, y)
        output_mask(x, y) = 2;
        % end
    end
end
close(h);
end
