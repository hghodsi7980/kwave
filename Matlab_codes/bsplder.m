function x = bsplder(x, k)

%     Y = BSPLDER(X, K) Returns the control points of the 'K'-th derivative
%     'Y' of the Bézier curve 'X'.

    % Calculate curve order
    n = size(x, 1);
    
    % K-differenciation
    x = diff(x, k) * prod(n - k : n - 1);
end