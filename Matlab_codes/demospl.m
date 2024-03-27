% Control points
x = rand(randi(2) + 2, 2);

% Open composite curve
y = bsplpcw(x);

% Closed composite curve
z = bsplpcw(x, []);

% Curve evaluations
a = bspl(x);
b = bspl(y);
c = bspl(z);

% Pre-allocate graphical objects
figure
tiledlayout(1, 3)
nexttile
hold on

% Plot Bézier curve
plot(a(:, 1), a(:, 2))
plot(x(:, 1), x(:, 2), 'o-')
pbaspect([1, 1, 1])
daspect([1, 1, 1])
bplt

% Plot open piece-wise curve
nexttile
hold on
plot(b(:, 1), b(:, 2))
plot(y(:, 1), y(:, 2), 'o-')
pbaspect([1, 1, 1])
daspect([1, 1, 1])
bplt

% Plot closed piece-wise curve
nexttile
hold on
plot(c(:, 1), c(:, 2))
plot(z(:, 1), z(:, 2), 'o-')
pbaspect([1, 1, 1])
daspect([1, 1, 1])
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