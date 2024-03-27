function y = bez2ber(x)

%     Y = BEZ2BER(X, DIM) Returns the Bernstein or polynomial coefficients
%     of the Bézier curve of control points 'X' along the dimension 'DIM'.

    % Parse inputs
    if nargin > 1 && ~isempty(dim) && dim == 1
        x = x';
    end
    
    % Initialise variables
    n = size(x, 1) - 1;
    b = cumprod(1 : n);
    
    % Lowest order point is shared
    y = x;
    
    % Loop over control points
    for i = 1 : n - 1
        s = (-1) ^ i;
        k = s * x(1, :) / b(i);
        
        % Recursive factorial division from binomial expansion
        for j = 1 : i - 1
            s = -s;
            k = k + s * x(j + 1, :) / (b(j) * b(i - j));
        end
        k = k + x(i + 1, :) / b(i);
        
        % Get Bernstein point of matching binomial order
        y(i + 1, :) = prod(n - i + 1 : n) * k;
    end
    
    % Get final point using the same approach
    s = (-1) ^ n;
    k = s * x(1, :) / b(n);
    for j = 1 : n - 1
        s = -s;
        k = k + s * x(j + 1, :) / (b(j) * b(n - j));
    end
    y(n + 1, :) = b(n) * k + x(n + 1, :);
end