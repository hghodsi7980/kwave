function angles = calculate_angles_for_nodes(vectors)
    % Normalize the vectors
    vectors = vectors ./ vecnorm(vectors, 2, 2);  % Normalize the vectors
    
    % Initialize angles array
    angles = zeros(size(vectors, 1), 1);
    
    % Calculate angles between consecutive vectors
    for i = 1:size(vectors, 1)
        vec1 = vectors(i, :);
        if i < size(vectors, 1)
            vec2 = vectors(i + 1, :);
        else
            vec2 = vectors(1, :);
        end
        
        % Calculate the angle between vec1 and vec2
        dot_product = dot(vec1, vec2);
        angle = acos(dot_product); % angle in radians
        
        % Convert to degrees
        angles(i) = rad2deg(angle);
    end
    
    % Display sum of angles to check if it is close to 360
    sum_of_angles = sum(angles);
    fprintf('Sum of angles: %f degrees\n', sum_of_angles);
    
    % Display angles
    disp('Angles in degrees:');
    disp(angles);
end

