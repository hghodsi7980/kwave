function [val, idx] = bspldis(x, n, b)

%     [VAL, IDX] = BSPLDIS(X, N, B) Returns the values 'VAL' and knots
%     'IDX' of the Bézier curve of control points 'X' divided into 'N'
%     equally-spaced intervals. A resolution factor 'B' can optionally
%     be specified for the routine BSPLARC.


    % Default: split in half
    if nargin < 2 || isempty(n)
        n = 2;
    end
    if nargin < 3 || isempty(b)
        b = 10 * n * size(x, 1);
    end
    
    % Initialise variables
    idx = zeros(n + 1, 1);
    idx(n + 1) = 1;
    
    % Get approximate arc-length
    a = bsplarc(x, b) / n;
    
    % Get resolution of the calculations
    b = [b; max(fix(b / n), 3)];
    b = fix(b(1) : (b(2) - b(1)) / (n - 2) : b(2));
    
    % Loop per partition: cumulative arc-length
    for i = 1 : n - 1
        t = idx(i) : (1 - idx(i)) / (b(i) - 1) : 1;
        k = bsplarc(x, t, -1) - a;
        [~, k] = min(abs(k));
        idx(i + 1) = t(k);
    end
    
    % Evaluate curve at splitting point
    val = bspl(x, idx);
end