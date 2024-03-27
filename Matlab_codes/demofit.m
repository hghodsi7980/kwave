% Number of points to fit
n = 5;

% Dummy curve to place points
X = bspldis(rand(4, 2), n - 1);

% Affine-invariant distribution
tic
[A, t1] = bsplget(X);
toc

% Arc-length distribution
tic
[B, t2] = bsplget(X, 0);
toc

% Uniform knot distribution
tic
[C, t3] = bsplget(X, 1);
toc

% Minimum length distribution
tic
[D, t4] = bsplget(X, 1e3);
toc

% Fitted curve with lower order
tic
[E, ~, t5] = bsplfit(X, n - 2, 1);
toc

% Spline evaluation for visualisation
a = bspl(A);
b = bspl(B);
c = bspl(C);
d = bspl(D);
e = bspl(E);

% Pre-allocate graphical object
figure
tiledlayout(1, 2)
nexttile
hold on

% Allocate curves that do not diverge
h = cell(5, 1);
if max(abs(a(:))) < 2
    h{1} = 'A-I';
    plot(a(:, 1), a(:, 2), 'Color', '#0072BD')
end
if max(abs(b(:))) < 2
    h{2} = 'A-L';
    plot(b(:, 1), b(:, 2), 'Color', '#D95319')
end
if max(abs(c(:))) < 2
    h{3} = 'U';
    plot(c(:, 1), c(:, 2), 'Color', '#EDB120')
end
if max(abs(d(:))) < 2
    h{4} = 'M-L';
    plot(d(:, 1), d(:, 2), 'Color', '#7E2F8E')
end
if max(abs(e(:))) < 2
    h{5} = 'FIT';
    plot(e(:, 1), e(:, 2), 'Color', '#77AC30')
end
plot(X(:, 1), X(:, 2), 'o', 'Color', 'k')

% Format plot
title('Curve Fitting')
h = h(~cellfun('isempty', h));
legend(h);
pbaspect([1, 1, 1])
bplt

% Knot visualisation
nexttile
hold on

% Allocate curves that do not diverge
h = cell(5, 1);
if max(abs(a(:))) < 2
    h{1} = 'A-I';
    plot(t3, t1, '.-', 'Color', '#0072BD')
end
if max(abs(b(:))) < 2
    h{2} = 'A-L';
    plot(t3, t2, '.-', 'Color', '#D95319')
end
if max(abs(c(:))) < 2
    h{3} = 'U';
    plot(t3, t3, '.-', 'Color', '#EDB120')
end
if max(abs(d(:))) < 2
    h{4} = 'M-L';
    plot(t3, t4, '.-', 'Color', '#7E2F8E')
end
if max(abs(e(:))) < 2
    h{5} = 'FIT';
    plot(t3, t5, '.-', 'Color', '#77AC30')
end

% Format plot
title('Knot Vector')
h = h(~cellfun('isempty', h));
h = legend(h);
pbaspect([1, 1, 1])
bplt
set(h, 'Location', 'NorthWest')


% Format plot
function bplt()

    % Format axes
    h = findall(gca, 'Type', 'axes');
    set(h, 'box', 'on', 'Color', 'w')
    set(h, 'TickDir', 'in')
    set(h, 'TickLabelInterpreter', 'LaTeX')

    % Format text
    h = findall(gca, 'Type', 'Text');
    set(h, 'Interpreter', 'LaTeX')
    set(h, 'FontSize', 2)

    % Format title
    h = get(gca, 'Title');
    set(h, 'FontSize', 12)

    % Format legend
    h = get(gca, 'Legend');
    set(h, 'Box', 'off')
    set(h, 'FontSize', 8)
    set(h, 'Location', 'Best')
    set(h, 'Interpreter', 'LaTeX')
end