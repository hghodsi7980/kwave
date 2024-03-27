% Generate random control points
X = rand(4, 3);
Y = rand(5, 3);
Z = rand(6, 3);

% Check dimensions
a = size(X, 2);
b = size(Y, 2);
c = size(Z, 2);
if a ~= b || a ~= c || b ~= c
    error('Dimensional incompatibility')
end

% Get Bézier evaluations for mathematical operators
x = bspl(X);
y = bspl(Y);
z = bspl(Z);

% 'X - Y + Z' and its residual
a = bsplfun(X, -Y, Z, '+');
r1 = x - y + z - bspl(a);
r1 = sqrt(mean(sum(r1 .* r1, 2)));

% 'cross(X, Y) * Z' and its residual
b = bsplfun(X, Y, 'x');
if size(X, 2) == 2
    b = bsplfun(b, Z(:, 1), '*');
else
    b = bsplfun(b, Z, '*');
end
if size(X, 2) == 2
    r2 = x(:, 1) .* y(:, 2) - x(:, 2) .* y(:, 1);
    r2 = r2 .* z - bspl(b);
else
    r2 = cross(x, y, 2) .* z - bspl(b);
end
r2 = sqrt(mean(sum(r2 .* r2, 2)));

% 'X * Y ^ 2 * Z ^ 3' and its residual
c = bsplfun(X, Y, Y, Z, Z, Z, '*');
r3 = x .* y .^ 2 .* z .^ 3 - bspl(c);
r3 = sqrt(mean(sum(r3 .* r3, 2)));

% Pre-allocate graphical objects
figure
tiledlayout(1, 3)

% Plot 'X - Y + Z'
nexttile
plot(bspl(a));
pbaspect([1, 1, 1])
title('$$X - Y + Z$$')
bplt

% Plot 'cross(X, Y) * Z'
nexttile
plot(bspl(b))
pbaspect([1, 1, 1])
title('$$(X \times Y) \cdot Z$$')
bplt

% Plot 'X * Y ^ 2 * Z ^ 3'
nexttile
plot(bspl(c))
pbaspect([1, 1, 1])
title('$$X \cdot Y ^ 2 \cdot Z ^ 3$$')
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