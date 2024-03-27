function [y, idx, val] = bsplcsp(x, y)

%     [Y, IDX, VAL] = BSPLCSP(X, Y) Returns the convexity in terms of cusps
%     and inflections of 'X'. Parsing a second input variable 'Y' specifies
%     the calculation of the inflection points of 'X' (non-enpty 'Y' value)
%     or the calculation of the cusps (empty or null 'Y').
%
%     The output parameters are:
%     
%         - 'Y'         Coordinates of the cusps and inflections
%         - 'IDX'       Knot values for the points in 'Y'
%         - 'VAL'       Curvature value at the points in 'Y'


    % Special case: inflections
    if nargin > 1 && ~isempty(y)
        idx = bsplint(bsplfun(bsplkin(bsplder(x, 2), 1), bsplder(x, 1), ...
              'cross'), 0);
        val = 0 * idx;
        y = bspl(x, idx);
        return
    end

    % Initialise variables
    p = 1e4;
    q = 1e2;
    
    % Default: local concavities
    val = bsplcrv(x, p);
    if nargin == 1
        
        % Cusps and inflections
        val = abs(val);
    end
    
    % Find local extrema
    idx = diff(val);
    idx = find(idx(1 : p - 2) .* idx(2 : p - 1) < 0) + 1;
    if isempty(idx)
        y = [];
        val = [];
        return
    end
    y = val(idx);
    
    % Loop over the number of ocurrences
    for i = 1 : length(idx)
        j = idx(i);
        
        % knot span refinement
        t = [j - 5, j + 5] / p;
        t(t < 0) = 0;
        t(t > 1) = 1;
        idx(i) = j / p;
        res = 1;
        iter = 0;
        
        % Loop until numerical precision is reached
        while res > eps
            
            % Evaluate curvature over refined span
            t = t(1) : (t(2) - t(1)) / (q - 1) : t(2);
            val = abs(bsplcrv(x, t));
            j = diff(val);
            j = find(j(1 : q - 2) .* j(2 : q - 1) < 0) + 1;
            if length(j) ~= 1
                break
            end
            
            % Update residuals
            res = abs(val(j) - y(i));
            y(i) = val(j);
            idx(i) = t(j);
            t = [t(j - 1), t(j + 1)];
            
            % Exit condition
            iter = iter + 1;
            if iter == 5
                break
            end
        end
    end
    
    % Update outputs
    if nargin == 1
        val = y;
        y = bspl(x, idx);
    else
        for i = 1 : length(idx)
            z = abs(bsplcrv(x, idx(i) * (0.99 : 0.01 : 1.01)));
            if z(2) < z(1) || z(2) < z(3)
                idx(i) = 0;
            end
        end
        val = y(idx > 0, :);
        idx = idx(idx > 0);
        y = bspl(x, idx);
    end
end