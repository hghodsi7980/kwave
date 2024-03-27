% Control points of the first curve to analyse
A = [0.7, 0.75; 0, 0.5; 0.2, 0.8; 0.1, 0.4; 1, 0; 0, 0.25];

% Get Bézier curve for visualisation
a = bspl(A, 1e3);

% Get curvature on same query knots
k = bsplcrv(A, 1e3);

% Get normal vector
[~, n] = bsplvec(A, 1e3);

% Superimpose scaled curvature to Bézier curve
k = a + 0.1 * range(a) .* n .* k / max(abs(k));

% Find cusps and inflections (concavities and convexities)
[b, idx, val] = bsplcsp(A);

% Pre-allocate graphical objects
figure
tiledlayout(1, 2)
nexttile
hold on

% Plot Bézier curve and its superimposed curvature
plot(a(:, 1), a(:, 2))
plot(k(:, 1), k(:, 2))
if ~isempty(b)
    
    % Loop over number of occurrences
    for i = 1 : length(val)
        
        % Flag for cusps (otherwise: inflections)
        csp = val(i) > abs(bsplcrv(A, idx(i) + 0.02 * (rand - 0.5)));
        if csp
            plot(b(i, 1), b(i, 2), 'o', 'Color', '#7E2F8E')
        else
            plot(b(i, 1), b(i, 2), 'o', 'Color', '#EDB120')
        end
    end
end

% Format plot
pbaspect([1, 1, 1])
daspect([1, 1, 1])
bplt


% Control points of the first curve to analyse
A = [0.5, 0.25; 0.5, -0.4; 1, 0.5; 0.3, 0.3; 0.7, 0];

% Get Bézier curve for visualisation
a = bspl(A, 1e3);

% Get curvature on same query knots
k = bsplcrv(A, 1e3);

% Get normal vector
[~, n] = bsplvec(A, 1e3);

% Superimpose scaled curvature to Bézier curve
k = a + 0.1 * range(a) .* n .* k / max(abs(k));

% Find cusps (concavities)
[b, idx, val] = bsplcsp(A, []);

nexttile
hold on

% Plot Bézier curve and its superimposed curvature
plot(a(:, 1), a(:, 2))
plot(k(:, 1), k(:, 2))
if ~isempty(b)
    plot(b(:, 1), b(:, 2), 'o', 'Color', '#7E2F8E')
end

% Format plot
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