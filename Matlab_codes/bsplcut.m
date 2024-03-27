function y = bsplcut(x, t)

%     Y = BSPLCUT(X, T) Returns the control points of the composite Bézier
%     curve 'Y' that results from splitting the curve of control points 'X'
%     at the query knots in 't' (by default, set as 0.5). If the parameter
%     't' is of the same dimensions as 'X', then the curve 'X' is splitted
%     at the closest point to the query values in 't' through projections.


    % Default: cut at the mid knot
    if nargin < 2 || isempty(t)
        t = 0.5;
    end
    [n, p] = size(t);
    
    % Cut knot is a 3D point
    if p > 1
        t = bsplpro(x, t);
        
    % Equal length segments
    elseif n == 1 && t(1) > 1
        n = t + 1;
        [~, t] = bspldis(x, t + 1);
        t = t(2 : n);
        
    % Query splitting knots
    else
        t = t(t > 0 & t < 1);
        n = length(t) + 1;
    end
    
    % Initialise control points
    m = n;
    [n, p] = size(x);
    y = zeros(m * n, p);
    
    % First segmentation
    z = bspldiv(x, n, t(1));
    p = 1 : n;
    y(p, :) = z(p, :);
    z(p, :) = [];
    
    % Loop over the remaining of the cutting knots
    for i = 2 : m - 1
        z = bspldiv(z, n, bsplpro(z, bspl(x, t(i))));
        y(n * (i - 1) + 1 : n * i, :) = z(p, :);
        z(p, :) = [];
    end
    
    % Add last segment to the final set of points
    y(n * (m - 1) + 1 : n * m, :) = z;
end

% Division algorithm: modified De-Casteljau algorithm
function y = bspldiv(x, n, u)
    m = 2 * n;
    v = 1 - u;
    y(m, :) = x(n, :);
    y(1, :) = x(1, :);
    for i = 1 : n - 1
        j = n - i;
        x = v * x(1 : j, :) + u * x(2 : j + 1, :);
        y([i + 1, m - i], :) = x([1, j], :);
    end
end