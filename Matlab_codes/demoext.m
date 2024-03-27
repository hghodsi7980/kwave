% Get control points of function to evaluate
X = rand(5, 2);

% Obtain its first two derivatives
Y = bsplder(X, 1);
Z = bsplder(Y, 1);

% Get zeros of first derivative
if size(X, 1) > 2
    [a, A] = bsplint(Y, [min(Y(:, 1)), 0; max(Y(:, 1)), 0]);
    [b, B] = bsplint(Y, [0, min(Y(:, 2)); 0, max(Y(:, 2))]);
else
    A = [];
    B = [];
end

% Get zeros of second derivative
if size(X, 1) > 3
    [c, C] = bsplint(Z, [min(Z(:, 1)), 0; max(Z(:, 1)), 0]);
    [d, D] = bsplint(Z, [0, min(Z(:, 2)); 0, max(Z(:, 2))]);
else
    C = [];
    D = [];
end

% Pre-allocate graphical objects
figure
tiledlayout(1, 3)

% Plot curve and its extrema
nexttile
hold on
x = bspl(X, 1e3);
plot(x(:, 1), x(:, 2))

% Plot extrema markers if they exist
if ~isempty(A)
    a = bspl(X, a(:, 1));
    plot(a(:, 1), a(:, 2), 'o', 'Color', '#D95319')
end
if ~isempty(B)
    b = bspl(X, b(:, 1));
    plot(b(:, 1), b(:, 2), 'o', 'Color', '#EDB120')
end
if ~isempty(C)
    c = bspl(X, c(:, 1));
    plot(c(:, 1), c(:, 2), 'o', 'Color', '#7E2F8E')
end
if ~isempty(D)
    d = bspl(X, d(:, 1));
    plot(d(:, 1), d(:, 2), 'o', 'Color', '#77AC30')
end

% Format plot
pbaspect([1, 1, 1])
daspect([1, 1, 1])
xx = get(gca, 'xlim');
yy = get(gca, 'ylim');
xx = xx + 0.1 * mean(range(x)) * [-1, 1];
yy = yy + 0.1 * mean(range(x)) * [-1, 1];
xlim(xx)
ylim(yy)
title('Extrema')
bplt


% Get first derivative curve
if ~isempty(Y)
    x = bspl(Y, 1e3);
else
    x = [0, 0];
end

% Plot first derivative and its zeros
nexttile
hold on
plot(x(:, 1), x(:, 2))
if ~isempty(A)
    plot(A(:, 1), A(:, 2), 'o', 'Color', '#D95319')
end
if ~isempty(B)
    b = bspl(Y, B(:, 1));
    plot(B(:, 1), B(:, 2), 'o', 'Color', '#EDB120')
end

% Format plot
pbaspect([1, 1, 1])
daspect([1, 1, 1])
xx = get(gca, 'xlim');
yy = get(gca, 'ylim');
xx = xx + 0.1 * mean(range(x)) * [-1, 1];
yy = yy + 0.1 * mean(range(x)) * [-1, 1];

% Plot lines at zeros
if ~isempty(A)
    p = plot(xx, [0, 0], 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', ':');
    p.Color(4) = 0.6;
end
if ~isempty(B)
    q = plot([0, 0], yy, 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', ':');
    q.Color(4) = 0.6;
end

% Format plot
xlim(xx)
ylim(yy)
title('$$d^1(X, Y) = 0$$')
bplt


% Get second derivative curve
if ~isempty(Z)
    x = bspl(Z, 1e3);
else
    x = [0, 0];
end

% Plot second derivative and its zeros
nexttile
hold on
plot(x(:, 1), x(:, 2))
if ~isempty(C)
    plot(C(:, 1), C(:, 2), 'o', 'Color', '#7E2F8E')
end
if ~isempty(D)
    plot(D(:, 1), D(:, 2), 'o', 'Color', '#77AC30')
end

% Format plot
pbaspect([1, 1, 1])
daspect([1, 1, 1])
xx = get(gca, 'xlim');
yy = get(gca, 'ylim');
xx = xx + 0.1 * mean(range(x)) * [-1, 1];
yy = yy + 0.1 * mean(range(x)) * [-1, 1];

% Plot lines at zeros
if ~isempty(C)
    p = plot(xx, [0, 0], 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', ':');
    p.Color(4) = 0.6;
end
if ~isempty(D)
    q = plot([0, 0], yy, 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', ':');
    q.Color(4) = 0.6;
end

% Format plot
xlim(xx)
ylim(yy)
title('$$d^2(X, Y) = 0$$')
bplt


% Format figure
function bplt()

    % Format axes
    h = findall(gca, 'Type', 'axes');
    set(h, 'box', 'on', 'Color', 'w')
    set(h, 'TickDir', 'in')
    set(h, 'TickLabelInterpreter', 'LaTeX')

    % Format text
    h = findall(gca, 'Type', 'Text');
    set(h, 'Interpreter', 'LaTeX')
    set(h, 'FontSize', 10)

    % Format title
    h = get(gca, 'Title');
    set(h, 'FontSize', 12)

    % Format legend
    h = get(gca, 'Legend');
    set(h, 'Box', 'off')
    set(h, 'Location', 'Best')
    set(h, 'Interpreter', 'LaTeX')
end