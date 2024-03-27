function y = bsplpcw(x, ~)

%     Y = BSPLPCW(X) Returns the composite control points of the set of
%     cubic Bézier curves 'Y' that interpolate the points in 'X' whilst
%     preserving C2 continuity. Retrieved and modified from:
%
%     Peter Lee (2009). Smooth Curve through a Set of 2D Points


    % Parse input arguments
    if nargin > 1
        x = [x; x(1, :)];
    end
    [m, p] = size(x);
    n = m - 1;
    q = n - 1;
    a(n, p) = 0; 
    d = a;
    
    % Linear interpolation case
    if n == 1
        y = [x(1, :); (2 * x(1, :) + x(2, :)) / 3; (4 * x(1, :) + 2 * ...
             x(2, :)) / 3 - x(1, :); x(2, :)];
        return
    end
    
    % Get curve end-points from interpolants
    t(n, :) = 4 * x(n, :) + 0.5 * x(m, :);
    t(1, :) = x(1, :) + 2 * x(2, :);
    t(2 : q, :) = 4 * x(2 : q, :) + 2 * x(3 : n, :);
    
    % Initialise cubic segments
    c = 2 * ones(1, p);
    a(1, :) = 0.5 * t(1, :);
    for i = 2 : q
        d(i, :) = 1 ./ c;
        c = 4 - d(i, :);
        a(i, :) = (t(i, :) - a(i - 1, :)) ./ c;
    end
    d(n, :) = 1 ./ c;
    c = 3.5 - d(n, :);
    
    % Calculate cubic segments
    a(n, :) = (t(n, :) - a(q, :)) ./ c;
    for i = 1 : q
        a(i, :) = a(i, :) - d(i + 1, :) .* a(i + 1, :);
    end
    b(n, :) = 0.5 * (x(m, :) + a(n, :));
    b(1 : q, :) = 2 * x(2 : n, :) - a(2 : n, :);
    
    % Build concatenated control point matrix
    ii = 1 : 4;
    y(4 * n - 3 : 4 * n, :) = [x(n, :); a(n, :); b(n, :); x(m, :)];
    for i = 1 : n - 1
        y(ii, :) = [x(i, :); a(i, :); b(i, :); x(i + 1, :)];
        ii = ii + 4;
    end
    
    % Correct closed case
    if norm(x(1, :) - x(m, :)) < 1e-6
        q = 4 * n - 1;
        a = y(2, :);
        b = y(q, :);
        c = 2 * x(1, :);
        while norm(c - a - b) > 1e-6
            b = 0.5 * (b + c - a);
            a = 0.5 * (a + c - b);
        end
        y(2, :) = a;
        y(q, :) = b;
    end
end