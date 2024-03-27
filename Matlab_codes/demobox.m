% Get control points for 2D sample
A = rand(5, 2);

% Get bézier curve for visualisation
B = bspl(A);

% Get tight bounding box
C = bsplbox(A, [], 1);

% Get standard bounding box
A = bsplbox(A, [], 0);

% Pre-allocate graphical objects
figure
tiledlayout(1, 2)
nexttile
hold on

% Plot curve
plot(B(:, 1), B(:, 2), 'Color', '#0072BD')

% Plot tight box's polygon
plot(C([1 : 4, 1], 1), C([1 : 4, 1], 2), 'Color', '#D95319')

% Plot standard box's polygon
plot(A([1 : 4, 1], 1), A([1 : 4, 1], 2), 'Color', '#EDB120')

% Format plot
pbaspect([1, 1, 1])
daspect([1, 1, 1])
title('Two-Dimensional')
bplt


% Get control points for 3D sample
D = rand(6, 3);
save = round(D, 2);

% Get Bézier curve for visualisation
E = bspl(D);

% Get tight bounding box
F = bsplbox(D, [], 1);

% Get standard bounding box
D = bsplbox(D, [], 0);

nexttile
hold on

% Plot curve
plot3(E(:, 1), E(:, 2), E(:, 3), 'Color', '#0072BD')

% Plot tight box's polygon
plot3(F([1, 2, 3, 4, 1], 1), F([1, 2, 3, 4, 1], 2), ...
      F([1, 2, 3, 4, 1], 3), 'Color', '#D95319')
plot3(F([5, 6, 7, 8, 5], 1), F([5, 6, 7, 8, 5], 2), ...
      F([5, 6, 7, 8, 5], 3), 'Color', '#D95319')
plot3(F([1, 2, 6, 5, 1], 1), F([1, 2, 6, 5, 1], 2), ...
      F([1, 2, 6, 5, 1], 3), 'Color', '#D95319')
plot3(F([4, 3, 7, 8, 4], 1), F([4, 3, 7, 8, 4], 2), ...
      F([4, 3, 7, 8, 4], 3), 'Color', '#D95319')
  
% Plot standard box's polygon
plot3(D([1, 2, 3, 4, 1], 1), D([1, 2, 3, 4, 1], 2), ...
      D([1, 2, 3, 4, 1], 3), 'Color', '#EDB120')
plot3(D([5, 6, 7, 8, 5], 1), D([5, 6, 7, 8, 5], 2), ...
      D([5, 6, 7, 8, 5], 3), 'Color', '#EDB120')
plot3(D([1, 2, 6, 5, 1], 1), D([1, 2, 6, 5, 1], 2), ...
      D([1, 2, 6, 5, 1], 3), 'Color', '#EDB120')
plot3(D([4, 3, 7, 8, 4], 1), D([4, 3, 7, 8, 4], 2), ...
      D([4, 3, 7, 8, 4], 3), 'Color', '#EDB120')

% Format plot
view(3)
pbaspect([1, 1, 1])
daspect([1, 1, 1])
title('Three-Dimensional')
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