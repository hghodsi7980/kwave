function k = bsplcrv(x, u)

%     K = BSPLCRV(X, U) calculates the curvature 'K' of a Bézier curve with
%     control points 'X' in the knot vector or knot resolution in 'U'.


    % Parse inputs
    if nargin < 2 || isempty(u)
        u = (0 : 1 / 999 : 1)';
    end
    [a, b] = size(x);

    % Get curve derivatives over query knot values
    k = bsplder(x, 2);
    if a == 3
        k = k .* ones(max(u(1), length(u)), 1);
    else
        k = bspl(k, u);
    end
    u = bspl(bsplder(x, 1), u);
    
    % Calculate cross product norm
    if b == 2
        k = u(:, 1) .* k(:, 2) - u(:, 2) .* k(:, 1);
    else
        k = cross3d(u, k);
        k = k .* k;
        k = sum(k, 2);
        k = sqrt(k);
    end
    
    % Obtain signed curvature
    u = u .* u;
    u = sum(u, 2);
    u = sqrt(u);
    u = u .* u .* u;
    k = k ./ u;
end

% Faster cross-product than cross(A, B, dim)
function c = cross3d(a, b)
    c(:, 3) = a(:, 1) .* b(:, 2) - a(:, 2) .* b(:, 1);
    c(:, 1) = a(:, 2) .* b(:, 3) - a(:, 3) .* b(:, 2);
    c(:, 2) = a(:, 3) .* b(:, 1) - a(:, 1) .* b(:, 3);
end