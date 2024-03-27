function [idx, val, res] = bsplpro(x, y)

%     [VAL, IDX, RES] = BSPLPRO(X, Y) Returns the values 'VAL', knots 'IDX'
%     and residuals 'RES' of the projections of the points 'Y' to the curve
%     generated from the control points in 'X'.


    % Initialise variables
    n = size(y, 1);
    p = 1e4;
    q = 1e2;
    idx = zeros(n, 1);
    val = bspl(x, p);
    
    % Loop over query points
    for i = 1 : n
        
        % Find closest knot
        j = val - y(i, :);
        j = j .* j;
        j = sum(j, 2);
        j = sqrt(j);
        [j, k] = min(j);
        t = k / p;
        
        % Get a target knot span
        k = [k - 5, k + 5] / p;
        k(k < 0) = 0;
        k(k > 1) = 1;
        res = 1;
        iter = 0;
        
        % Iterate until numerical precision is met
        while res > 1e-6
            
            % Evaluate minimum distance
            u = k(1) : (k(2) - k(1)) / (q - 1) : k(2);
            jj = bspl(x, u) - y(i, :);
            jj = jj .* jj;
            jj = sum(jj, 2);
            jj = sqrt(jj);
            [jj, kk] = min(jj);
            
            % Decreasing residual: get knot span
            if jj < j
                res = abs(jj - j);
                j = jj;
                t = u(kk);
                kk = [kk - 1, kk + 1];
                kk(kk < 0) = 0;
                kk(kk > 1) = 1;
                k = u(kk);
                q = 1e2;
                
            % Otherwise: increase knot resolution
            else
                q = q + 1;
            end
            
            % Exit condition
            iter = iter + 1;
            if iter == 1e2
                break
            end
        end
        idx(i) = t;
    end
    
    % Evaluate goodness of projection
    val = bspl(x, idx);
    if nargout > 2
        res = val - y;
        res = res .* res;
        res = sum(res, 2);
        res = sqrt(res);
    end
end