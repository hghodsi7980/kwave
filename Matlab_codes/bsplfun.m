function x = bsplfun(varargin)

%     X = BSPLFUN(X1, X2, ..., XN, 'OPER.') Applies the operator 'OPER.' to
%     the Bézier curves in 'X1', 'X2', ..., 'XN', and returns the resulting
%     Bézier curve in the form of its control points 'X'. All curves should
%     have the same dimension (number of columns 'q').
%
%     Examples of use are bsplfun(X, Y, '*') or bsplfun(X, -Y, Z, '+')


    % Initialise inputs
    n = numel(varargin) - 1;
    p = size(varargin{1}, 2);
    q = cellfun(@numel, varargin);
    q = q(1 : n) / p;
    fun = varargin{n + 1};
    
    % Summation of curves
    if strcmp(fun, '+') || strcmp(fun, 'plus')
        p = max(q);
        x = bsplkin(varargin{1}, p - q(1));
        for i = 2 : n
            x = x + bsplkin(varargin{i}, p - q(i));
        end
        return
    end
    
    % Product of curves
    if strcmp(fun, '*') || strcmp(fun, 'times')
        k = cumsum(q) - (0 : n - 1);
        b = binomial(k(n))';
        x = zeros(k(n), p);
        x(1 : k(1), :) = varargin{1} .* b(1 : k(1), k(1));
        for i = 2 : n
            y = varargin{i} .* b(1 : q(i), q(i));
            for j = p : -1 : 1
                x(1 : k(i), j) = conv(x(1 : k(i - 1), j), y(:, j));
            end
        end
        x = x ./ b(:, k(n));
        return
    end
    
    % Cross-product of curves
    if strcmp(fun, 'x') || strcmp(fun, 'cross')
        if p == 2 && n == 2
            x = bsplfun(varargin{1}(:, 1), varargin{2}(:, 2), '*') - ...
                bsplfun(varargin{1}(:, 2), varargin{2}(:, 1), '*');
        else
            k = cumsum(q) - (0 : n - 1);
            x = zeros(k(n), 3);
            x(1 : k(1), 1 : p) = varargin{1};
            for i = 2 : n
                y = [varargin{i}, zeros(3 - p)];
                j = 1 : k(i - 1);
                x(1 : k(i), :) = [bsplfun(x(j, 2), y(:, 3), '*') - ...
                                  bsplfun(x(j, 3), y(:, 2), '*'), ...
                                  bsplfun(x(j, 3), y(:, 1), '*') - ...
                                  bsplfun(x(j, 1), y(:, 3), '*'), ...
                                  bsplfun(x(j, 1), y(:, 2), '*') - ...
                                  bsplfun(x(j, 2), y(:, 1), '*')];
            end
            x = x(:, 1 : p);
        end
        return
    end
end

% Stable construction of the binomial terms (Pascal Triangle)
function b = binomial(n)
    b = zeros(n);
    for i = 1 : n
        j = i - 1;
        b(i, 1) = 1;
        b(i, i) = 1;
        for k = 2 : j
            b(i, k) = b(j, k - 1) + b(j, k);
        end
    end
end