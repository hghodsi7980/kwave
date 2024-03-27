function y = bsplarc(x, t, s)

%     L = BSPLARC(X, T, S) Returns the arc length 'L' of a Bézier curve of
%     control points 'X' using the settings in 's'. Syntax:
% 
%     L = BSPLARC(X)            Approximate total length with 1e6 points
%     L = BSPLARC(X, N)         Approximate total length with 'N' points
%     L = BSPLARC(X, [a, b])    Apporximate [a b] length with 1e6 points
%     L = BSPLARC(X, T, 1)      Analytic length given 'T'
%     L = BSPLARC(X, T, -1)     Approximate cumulative length given 'T'


    % Parse inputs
    if nargin < 2 || isempty(t)
        t = [0, 1];
    end
    if nargin < 3 || isempty(s)
        s = 0;
    end
    
    % Symbolic: exact derivation in Bernstein form
    if s == 1
        x = bez2ber(bsplder(x, 1));
        syms u
        y = 0 * x(1, :);
        for i = 1 : size(x, 1)
            y = y + x(i, :) * u ^ (i - 1);
        end
        y = double(int(norm(y), t));
        
    % Approximative results via discrete segments
    else
        if length(t) == 2
            y = 1e6;
            t = t(1) : (t(2) - t(1)) / (y - 1) : t(2);
        end
        y = bspl(x, t);
        y = diff(y);
        y = y .* y;
        y = sum(y, 2);
        y = sqrt(y);
        
        % Arc-length distribution or total arc-length
        if s == -1
            y = [0; cumsum(y)];
        else
            y = sum(y);
        end
    end
end