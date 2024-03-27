function y = bspldcj(x, u)

%     Y = BSPLDCJ(X, U) Returns the Bézier curve evaluation 'Y' from a set
%     of control points 'X' and a knot vector or resolution 'U' that is by
%     default set as 100 equally-spaced knots. The function is the simpler
%     version of the Bézier curve algorithm, and is affine to BSPL.


    % Parse input arguments
    [n, dim] = size(x);
    if nargin < 2 || isempty(u)
        u = (0 : 1 / 99 : 1)';
        m = 100;
    else
        u = double(u(:));
        m = length(u);
        if m == 1
            if u == 0
                y = x(1, :);
                return
            elseif u == 1
                y = x(n, :);
                return
            end
            if u == fix(u)
                m = u;
                u = (0 : 1 / (u - 1) : 1)';
            end
        else
            if any(u < 0)  || any(u > 1)
                u = u - min(u);
                u = u / max(u);
            end
        end
    end
    
    % Loop over the query knots
    v = 1 - u;
    y(m, dim) = 0;
    for i = 1 : m
        
        % Recursive interpolation
        t = x;
        for j = 1 : n
           for k = 1 : n - j
               t(k, :) = v(i) * t(k, :) + u(i) * t(k + 1, :);
           end
        end
        y(i, :) = t(1, :);
    end
end