function x = bsplkin(x, k)

    % Y = BSPLKIN(X, K) Returns the control points 'Y' of the Bézier curve
    % that results of inserting 'K' knots to the control points 'X'.

    
    % Initialise variables
    n = size(x, 1);
    a = x(1, :);
    b = x(n, :);
    n = n - 1;
    
    % Recursive addition of terms following de boor
    for n = n : n + k - 1
        x = [a; (x(1 : n, :) .* (1 : n)' + x(2 : n + 1, :) .* ...
            (n : -1 : 1)') / (n + 1); b];
    end
end