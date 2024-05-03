function volume = tetrahedron_volume(A, B, C, D)
    % Calculate vectors
    AB = B - A;
    AC = C - A;
    AD = D - A;

    % Calculate scalar triple product
    triple_product = dot(AB, cross(AC, AD));

    % Calculate volume
    volume = abs(triple_product) / 6;
end