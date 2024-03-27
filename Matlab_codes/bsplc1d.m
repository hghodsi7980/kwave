function [c, n] = bsplc1d(x, ~)

%     [C, N] = BSPLC1D(X, ~) Returns the discrete 1D curvature 'K' and the
%     normals 'N' of a discrete function 'X'. A second input specifies the
%     calculation of the signed curvature.


    % Parse input curve
    [n, p] = size(x);
    
    % Get finitie differences
    u = (x(2, :) - x(1, :)) ./ (x(3, :) - x(2, :));
    v = (x(n, :) - x(n - 1, :)) ./ (x(n - 1, :) - x(n - 2, :));
    
    % Prepare cross-products
    x = [x, zeros(n, 3 - p)];
    i = x(1 : n - 2, :) - x(3 : n, :);
    j = x(2 : n - 1, :) - x(3 : n, :);
    
    % Pre-calculate cross-products
    a = cross3d(i, j);
    b = cross3d(j, a);
    x = cross3d(i, a);
    
    % Initialise circumcentre terms
    i = i .* i;
    j = j .* j;
    a = a .* a;
    a = 2 * sum(a, 2);
    k = sum(i, 2) .* b;
    k = k - sum(j, 2) .* x;
    k = k ./ a;
    k = k(:, 1 : p);
    
    % Calculate circumcentre
    a = u .* (k(2, :) - k(1, :));
    b = v .* (k(n - 2, :) - k(n - 3, :));
    k = [k(1, :) - a; k; k(n - 2, :) + b];
    
    % Signed curvature
    if nargin > 1
        i = k(2 : n, :) .* k(1 : n - 1, :) < 0;
        i = all(i, 2);
        i = find(i) + 1;
        i = [1; i; n];
        
        % Get discrete curvature values
        c = k .* k;
        c = sum(c, 2);
        c = sqrt(c);
        c(isnan(c)) = Inf;
        c = 1 ./ c;
        
        % Correct sign
        for j = 2 : size(i, 1)
            c(i(j - 1) : i(j)) = c(i(j - 1) : i(j)) * (-1)^ j;
        end
        if -min(c) > max(c)
            c = -c;
        end
    else
        
        % Discrete curvature values
        c = k .* k;
        c = sum(c, 2);
        c = sqrt(c);
        c(isnan(c)) = Inf;
        c = 1 ./ c;
    end
    
    % Normals to curve
    if nargout > 1
        n = c .* k;
    end
end

% Faster cross-product than cross(A, B, dim)
function c = cross3d(a, b)
    c(:, 3) = a(:, 1) .* b(:, 2) - a(:, 2) .* b(:, 1);
    c(:, 2) = a(:, 3) .* b(:, 1) - a(:, 1) .* b(:, 3);
    c(:, 1) = a(:, 2) .* b(:, 3) - a(:, 3) .* b(:, 2);
end