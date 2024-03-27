function [u, v] = bsplvec(x, u, col)

%     [T, N] = BSPLVEC(X, U, COL) Calculates the tangent 'T' and normal 'N'
%     vectors to the Bézier curve of control points 'X' over the knots 'U',
%     set by default to 100 uniform points between 0 and 1. If 'COL' is set
%     to '1', the tangent and normal vectors are calculated for a matrix of
%     one-dimensional Bézier curves.


    % Parse input arguments
    if nargin < 2 || isempty(u)
        u = 100;
    end
    if nargin < 3 || isempty(col)
        col = 0;
    end
    [n, p] = size(x);
    if p == 1 || p > 3
        col = 1;
    end
    
    % Single curve
    if col == 0
        
        % Two-dimensional curves
        if p == 2
            u = bspl(bsplder(x, 1), u)';
            u = (u ./ sqrt(sum(u .* u)))';
            v = [-u(:, 2), u(:, 1)];
        else
            
            % Only tangent vector required
            if nargout < 2
                u = bspl(bsplder(x, 1), u)';
                u = (u ./ sqrt(sum(u .* u)))';
                
            % Build 3D tangents and normals
            else
                v = bspl(bsplder(x, 2), u)';
                u = bspl(bsplder(x, 1), u)';
                u = u ./ sqrt(sum(u .* u));
                v = cross(u + v, u);
                v = cross(v, u);
                v = v ./ sqrt(sum(v .* v));
                u = u';
                v = v';
            end
        end
        
    % Array of 1D curves: get tangents and normals in matrix form
    else
        x = reshape([repmat((0 : 1 / (n - 1) : 1), p, 1)'; x], n, []);
        u = bspl(bsplder(x, 1), u)';
        u = (u ./ sqrt(sum(u .* u)))';
        v = u(:, reshape([2 : 2 : 2 * p; 1 : 2 : 2 * p - 1], 1, []));
        v = -v .* (2 * rem(1 : 2 * p, 2) - 1);
    end
end