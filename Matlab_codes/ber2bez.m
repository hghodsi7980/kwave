function y = ber2bez(x, dim)

%     Y = BER2BEZ(X, DIM) Returns the Bézier coefficients 'Y' of the curve
%     in Bernstein or polynomial form with coefficients 'X'. The parameter
%     'DIM' optionally specifies the dimension to operate along.


    % Parse input arguments
    if nargin > 1 && ~isempty(dim) && dim == 1
        x = x';
    end
    
    % Initialise variables
    a = size(x, 1);
    b = ones(a, 1);
    c = a + 1;
    n = a - 1;
    
    % Obtain binomial form of the control points
    for i = 2 : n
        if i > c / 2
            b(i) = b(c - i);
        else
            b(i) = n;
            for j = 2 : i - 1
               b(i) = b(i) * (a - j) / j;
            end
        end
    end
    x = x ./ b;
    y = x;
    
    % Get Bézier points from forward sums
    for i = 1 : n
        x = x(1 : a - i, :) + x(2 : c - i, :);
        y(i + 1, :) = x(1, :);
    end
end