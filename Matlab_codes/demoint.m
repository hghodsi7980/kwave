% Control points of intersecting curves
A = [0, 0.2; 1.85, 0.85; -0.85, 0.85; 1, 0.2];
B = [0, 0.8; 1.85, 0.15; -0.85, 0.15; 1, 0.8];

% Bézier curves for visualisation
a = bspl(A, 1e3);
b = bspl(B, 1e3);

% Pre-allocate graphical objects
figure
tiledlayout(1, 3)
nexttile
hold on

% Plot Bézier curves
plot(a(:, 1), a(:, 2))
plot(b(:, 1), b(:, 2))

% Evaluate mutual intersections
[~, val] = bsplint(A, B);
if ~isempty(val)
    plot(val(:, 1), val(:, 2), 'o', 'Color', '#EDB120')
end

% Evaluate self-intersections
[~, val] = bsplint(A);
if ~isempty(val)
    plot(val(:, 1), val(:, 2), 'o', 'Color', '#EDB120')
end
[~, val] = bsplint(B);
if ~isempty(val)
    plot(val(:, 1), val(:, 2), 'o', 'Color', '#EDB120')
end

% Format plot
ylim([0.2, 0.8])
pbaspect([1 1 1])
daspect([1 1 1])
bplt


% Control points of second set of curves
A = [0.25, -1; 0.4, 2.5; 0.6, -2.5; 0.75, 1];
B = 0.5 - A(4 : -1 : 1, 2 : -1 : 1);

% Bézier curves for visualisation
a = bspl(A, 1e3);
b = bspl(B, 1e3);

nexttile
hold on

% Plot Bézier curves
plot(a(:, 1), a(:, 2))
plot(b(:, 1), b(:, 2))

% Evaluate mutual intersections
[~, val] = bsplint(A, B);
if ~isempty(val)
    plot(val(:, 1), val(:, 2), 'o', 'Color', '#EDB120')
end

% Evaluate self-intersections
[~, val] = bsplint(A);
if ~isempty(val)
    plot(val(:, 1), val(:, 2), 'o', 'Color', '#EDB120')
end
[~, val] = bsplint(B);
if ~isempty(val)
    plot(val(:, 1), val(:, 2), 'o', 'Color', '#EDB120')
end

% Format plot
xlim([0, 1])
pbaspect([1 1 1])
daspect([1 1 1])
bplt



% Control points of the third demo
A = [0, 0; 0.6, 0.45; -4.4, 1.55; 5.65, -12.15; 4.45, 21; ...
    -12.25, -17.35; 12.35, 12.8; -12.25, -5.15; 4.65, -1.85; ...
    -0.6, 0.1; 0.95, -1.20; 1, 0];

nexttile
hold on

% Plot Bézier curve
a = bspl(A, 1e3);
plot(a(:, 1), a(:, 2))
x = xlim;
y = ylim;

% Evaluate intersections with (x, y) = 0
[~, val] = bsplint(A, 0);
if ~isempty(val)
    if any(round(val(:, 1), 2) == 0)
        plot([0, 0], y, 'Color', '#D95319')
    end
    if any(round(val(:, 2), 2) == 0)
        plot(x, [0, 0], 'Color', '#D95319')
    end
    plot(val(:, 1), val(:, 2), 'o', 'Color', '#EDB120')
end

% Format plot
pbaspect([1 1 1])
xlim(x)
ylim(y)
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