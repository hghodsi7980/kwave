function curve = bezier_volume(control_points, t, s, u)
    % Control points should be a 3D array [n+1, m+1, p+1, 3]
    % where n, m, and p are the degrees of the curve in each dimension
    % t, s, and u are the parameter values in each dimension

    % Initialize the volume
    curve = zeros(numel(t), numel(s), numel(u), 3);

    % Get the degrees
    n = size(control_points, 1) - 1;
    m = size(control_points, 2) - 1;
    p = size(control_points, 3) - 1;

    % Compute basis functions
    Bt = bernstein_basis(n, t);
    Bs = bernstein_basis(m, s);
    Bu = bernstein_basis(p, u);

    % Compute the volume using tensor product of basis functions and control points
    for i = 1:numel(t)
        for j = 1:numel(s)
            for k = 1:numel(u)
                for ii = 1:(n+1)
                    for jj = 1:(m+1)
                        for kk = 1:(p+1)
                            curve(i, j, k, :) = curve(i, j, k, :) + ...
                                Bt(i, ii) * Bs(j, jj) * Bu(k, kk) * control_points(ii, jj, kk, :);
                        end
                    end
                end
            end
        end
    end
end

function B = bernstein_basis(n, t)
    % Compute Bernstein basis functions of degree n at parameter values t
    % Returns an array of size [numel(t), n+1]

    % Initialize the array for basis functions
    B = zeros(numel(t), n + 1);

    % Compute basis functions using recursive formula
    for i = 0:n
        B(:, i + 1) = nchoosek(n, i) * t.^i .* (1 - t).^(n - i);
    end
end
