function y = bsplbox(x, dim, tight)

%     Y = BSPLBOX(X, DIM, TIGHT) Generates a bounding box around the Bézier
%     curve of control points 'X' in the dimension 'DIM' (by default, 2). A
%     third optional parameter 'TIGHT' queries the calculation of the tight
%     (default) or the standard bounding box. The matrix 'Y' has the vertex
%     coordinates of the bounding box.


    % Parse inputs
    if nargin < 2 || isempty(dim)
        dim = 2;
    end
    if nargin < 3 || isempty(tight)
        tight = 1;
    end
    p = size(x, dim);
    
    % Special cases: single or high-dimensional orders
    if p == 1 || p > 3
        a = 1e4;
        x = bspl(x, a);
        y = [min(x); max(x)]';
    elseif p == 2
        
        % Bounding box
        if tight == 0
            a = 1e4;
            x = bspl(x, a);
            y = [min(x); max(x)];
            y = reshape(y([1, 2, 2, 1, 3, 3, 4, 4]), [], 2);
            return
        end
        
        % Tight bounding box
        a = 1e2;
        b = 1e2;
        t = 2 * pi * (0 : 1 / (b - 1) : 1);
        c1 = cos(t);
        s1 = sin(t);
        
        % Loop over query angles
        u = Inf;
        for i = 1 : b
            r = [c1(i), -s1(i); s1(i), c1(i)];
            z = bspl(x, a) * r;
            z = [min(z); max(z)];
            s = prod(diff(z));
            if s < u
                u = s;
                y = z;
                t = r;
            end
        end
        
        % Reshape to vertex form
        y = reshape(y([1, 2, 2, 1, 3, 3, 4, 4]), [], 2) / t;
    else
        
        % Bounding box
        if tight == 0
            a = 1e4;
            x = bspl(x, a);
            y = [min(x); max(x)];
            y = reshape(y([1, 2, 2, 1, 1, 2, 2, 1, 3, 3, 4, 4, 3, 3, 4, ...
                4, 5, 5, 5, 5, 6, 6, 6, 6]), [], 3);
            return
        end
        
        % Tight bounding box
        a = 1e2;
        b = 10;
        m = b * b * b;
        t = zeros(m, 3);
        s = 2 * pi * (0 : 1 / (b - 1) : 1)';
        
        % Get rotation angles
        for i = 1 : 3
            r = [b, b, b];
            r(i) = 1;
            t(:, i) = reshape(repmat(permute(s, [2 : i, 1, i + 1 : 3]), ...
                      r), m, []);
        end
        
        % Pre-calculate rotation matrix terms
        c1 = cos(t(:, 1));
        c2 = cos(t(:, 2));
        c3 = cos(t(:, 3));
        s1 = sin(t(:, 1));
        s2 = sin(t(:, 2));
        s3 = sin(t(:, 3));
        
        % Construction of the rotation matrix given the query angles
        d1 = c2 .* c3;
        d2 = c2 .* s3;
        d3 = -s2;
        d4 = s1 .* s2 .* c3 - c1 .* s3;
        d5 = s1 .* s2 .* s3 + c1 .* c3;
        d6 = s1 .* c2;
        d7 = c1 .* s2 .* c3 + s1 .* s3;
        d8 = c1 .* s2 .* s3 - s1 .* c3;
        d9 = c1 .* c2;
        
        % Loop on the query angles
        v = Inf;
        for i = 1 : m
            r = [d1(i), d4(i), d7(i);
                 d2(i), d5(i), d8(i);
                 d3(i), d6(i), d9(i)];
            z = bspl(x, a) * r;
            z = [min(z); max(z)];
            s = prod(diff(z));
            if s < v
                v = s;
                y = z;
                t = r;
            end
        end
        
        % Reshape to vertex form
        y = reshape(y([1, 2, 2, 1, 1, 2, 2, 1, 3, 3, 4, 4, 3, 3, 4, 4, ...
            5, 5, 5, 5, 6, 6, 6, 6]), [], 3) / t;
    end
end